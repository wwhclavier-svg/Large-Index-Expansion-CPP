(* ::Package:: *)
(* VerifyExpand-Prepare.wl — IBP family definition + region solving *)
(* Usage: wolframscript -file VerifyExpand-Prepare.wl <famname> *)
(*
    Pipeline:
      Prepare → MMAExpand → Compare

    This file also serves as a shared utility module: later steps can
      Get["VerifyExpand-Prepare.wl"]
    to access paths, CLI helpers, file search, and formatting functions
    without re-executing the Prepare computation.
*)

SetDirectory[DirectoryName[$InputFileName]];
Get["LIEWorkflow.wl"];

(* ============================================================ *)
(* SECTION: Shared Utilities                                    *)
(* ============================================================ *)

(* ---- Paths ---- *)
$CPPVerifyRoot = ParentDirectory[$InputFileName] <> "/";
$MMARoot = DirectoryName[$InputFileName];

(* ---- CLI ---- *)
parseFamilyArg[] := Module[{famname},
    If[Length[$ScriptCommandLine] < 2,
        Print["ERROR: Missing <famname> argument."];
        Exit[1]
    ];
    $ScriptCommandLine[[2]]
];

(* ---- Family config from database ---- *)
loadFamilyConfig[famname_] := Module[{db},
    db = $CPPVerifyRoot <> "FamilyDatabase/FamilyDatabase.wl";
    If[!FileExistsQ[db],
        Print["ERROR: FamilyDatabase not found: ", db];
        Exit[1]
    ];
    Get[db];
    If[!KeyExistsQ[$FamilyDatabase, famname],
        Print["ERROR: Unknown family '", famname, "'."];
        Do[Print["  ", k], {k, Keys[$FamilyDatabase]}];
        Exit[1]
    ];
    $FamilyDatabase[famname]
];

(* ---- Target directory ---- *)
targetDir[famname_] := Module[{d},
    d = $CPPVerifyRoot <> famname <> "/";
    If[!FileExistsQ[d], CreateDirectory[d, CreateIntermediateDirectories -> True]];
    d
];

(* ---- File search across known locations ---- *)
findResultFile[famname_, baseName_] := Module[{paths, found},
    paths = {
        $CPPVerifyRoot <> famname <> "/",
        "./",
        "../../",
        "../../build/"
    };
    found = SelectFirst[paths, FileExistsQ[# <> baseName] &];
    If[MissingQ[found], baseName, found <> baseName]
];

(* ---- Formatting: compress polynomial for table display ---- *)
compressPoly[poly_, maxTerms_:5] := Module[{terms, nTerms},
    If[poly === 0, Return["0"]];
    terms = If[Head[poly] === Plus, List @@ poly, {poly}];
    nTerms = Length[terms];
    If[nTerms <= maxTerms,
        ToString[poly, InputForm],
        ToString[Plus @@ Take[terms, maxTerms], InputForm] <> " + ... (" <> ToString[nTerms] <> " terms)"
    ]
];

(* ---- Formatting: solution rules for Markdown ---- *)
formatSolution[rule_] := If[Length[rule] === 0,
    "(none)",
    StringRiffle[rule /. Rule[a_, b_] :> ToString[a] <> " -> " <> ToString[b, InputForm], ", "]
];

(* ============================================================ *)
(* SECTION: Main — Prepare                                      *)
(* ============================================================ *)

runPrepare[famname_] := Module[{
    config, data, expregdata, timeDefine, timeRegions,
    target, regionSummary, sectors, Alist, ne, char,
    binFile, ringFile, mmaTiming, checkpointFile
},
    config = loadFamilyConfig[famname];

    Print["========================================"];
    Print["VerifyExpand-Prepare     fam = ", famname];
    Print["  Description: ", config["Description"]];
    Print["  Numeric: ", config["Numeric"]];
    Print["  Modulus: ", config["Modulus"]];
    Print["========================================"];

    (* Step 1: Define Family *)
    Print["\n[Step 1] Defining IBP family..."];
    timeDefine = AbsoluteTiming[
        data = LIEDefineFamily[
            config["Propagators"], config["LoopMomenta"],
            config["ExternalMomenta"], config["KinematicRules"],
            config["TopSector"],
            "Numeric" -> config["Numeric"],
            Modulus -> config["Modulus"],
            Verbose -> False
        ];
    ][[1]];
    Print["  NE = ", data["Family", "NE"]];
    Print["  NIBP = ", data["Family", "NIBP"]];
    Print["  Time (Family Definition): ", timeDefine, " s"];

    (* Step 2: Solve Regions *)
    Print["\n[Step 2] Solving regions..."];
    timeRegions = AbsoluteTiming[
        data = LIESolveRegions[data, Verbose -> False];
    ][[1]];
    expregdata = data["Regions"];
    Print["  Sectors with regions: ", Length[Keys[expregdata]]];
    Print["  Time (Region Solving): ", timeRegions, " s"];

    target = targetDir[famname];

    (* Export Region Summary *)
    Print["\n[Step 3] Exporting region summary..."];
    regionSummary = {};
    sectors = Sort[Keys[expregdata]];
    Do[
        nRegs = Length[expregdata[sector]];
        Do[
            coordRing = expregdata[[Key[sector], i, Key["CoordinateRing"]]];
            varDeg = coordRing[[Key["VarDeg"]]];
            AppendTo[regionSummary, <|
                "Sector" -> sector,
                "RegionIndex" -> i,
                "VarDeg" -> varDeg,
                "NumSolutions" -> Times @@ varDeg,
                "VarIndep" -> coordRing[[Key["VarIndep"]]],
                "VarDep" -> coordRing[[Key["VarDep"]]],
                "VarRule" -> coordRing[[Key["VarRule"]]],
                "MinPoly" -> coordRing[[Key["MinPoly"]]],
                "Nb" -> Length[coordRing[[Key["MonomialBasis"]]]]
            |>];
        , {i, nRegs}];
    , {sector, sectors}];
    Put[regionSummary, target <> "Compare-RegionInfo-" <> famname <> ".m"];
    Print["  Exported ", Length[regionSummary], " region(s)"];

    (* Export binary matrices *)
    Print["\n[Step 4] Exporting binary matrices..."];
    Alist = data["Family", "AList"];
    ne = Length[Alist];
    char = config["Modulus"];
    binFile = target <> "IBPMat_" <> famname <> ".bin";
    ringFile = target <> "RingData_" <> famname <> ".bin";
    Get["ExportBinary_IBPMatrix.wl"];
    ExportBinaryIBPMatrix[binFile, expregdata, char];
    ExportBinaryRingData[ringFile, expregdata, Alist, ne, char];
    Print["  ", binFile];
    Print["  ", ringFile];

    (* Export timing *)
    mmaTiming = <|"FamilyDefinition" -> timeDefine, "RegionSolving" -> timeRegions|>;
    Export[target <> "Compare-MMATiming-" <> famname <> ".m", mmaTiming, "Text"];

    (* Save checkpoint for MMAExpand *)
    Print["\n[Step 5] Saving checkpoint..."];
    checkpointFile = target <> "PrepareCheckpoint-" <> famname <> ".wdx";
    Export[checkpointFile, data, "WDX"];
    Print["  ", checkpointFile];

    Print["\n========================================"];
    Print["VerifyExpand-Prepare complete!"];
    Print["  Next: wolframscript -file VerifyExpand-MMAExpand.wl ", famname];
    Print["========================================"];
];

(* ============================================================ *)
(* Entry point: only execute when invoked as the primary script *)
(* ============================================================ *)
If[
    Length[$ScriptCommandLine] >= 2 &&
    $ScriptCommandLine[[1]] =!= "wolframscript" &&
    StringContainsQ[$ScriptCommandLine[[1]], "VerifyExpand-Prepare.wl"],
    runPrepare[parseFamilyArg[]]
];
