(* ::Package:: *)
(* Compare-FamilyGenerate.wl — Generic IBP family .bin generator *)
(* Usage: wolframscript -file Compare-FamilyGenerate.wl <famname> *)
(* Reads config from verify/FamilyDatabase/FamilyDatabase.wl *)

SetDirectory[DirectoryName[$InputFileName]];
Get["LIEWorkflow.wl"];
Get["ExportBinary_IBPMatrix.wl"];
Get[ParentDirectory[DirectoryName[$InputFileName]] <> "/FamilyDatabase/FamilyDatabase.wl"];

(* Parse CLI argument *)
If[Length[$ScriptCommandLine] < 2,
    Print["ERROR: Missing family name argument."];
    Print["Usage: wolframscript -file Compare-FamilyGenerate.wl <famname>"];
    PrintAllFamilies[];
    Exit[1]
];
famname = $ScriptCommandLine[[2]];

(* Check if family exists *)
If[!KeyExistsQ[$FamilyDatabase, famname],
    Print["ERROR: Unknown family '", famname, "'."];
    PrintAllFamilies[];
    Exit[1]
];

config = $FamilyDatabase[famname];

Print["========================================"];
Print["Compare-FamilyGenerate     fam = ", famname];
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

(* Step 2b: Export Region Summary *)
Print["\n[Step 2b] Exporting region summary..."];
regionSummary = {};
sectors = Sort[Keys[expregdata]];
Do[
  nRegs = Length[expregdata[sector]];
  Do[
    coordRing = expregdata[[Key[sector], i, Key["CoordinateRing"]]];
    varDeg = coordRing[[Key["VarDeg"]]];
    entry = <| "Sector" -> sector,
       "RegionIndex" -> i,
       "VarDeg" -> varDeg,
       "NumSolutions" -> Times @@ varDeg,
       "VarIndep" -> coordRing[[Key["VarIndep"]]],
       "VarDep" -> coordRing[[Key["VarDep"]]],
       "VarRule" -> coordRing[[Key["VarRule"]]],
       "MinPoly" -> coordRing[[Key["MinPoly"]]],
       "Nb" -> Length[coordRing[[Key["MonomialBasis"]]]]
    |>;
    AppendTo[regionSummary, entry];
  , {i, nRegs}];
, {sector, sectors}];
targetDir = "../" <> famname <> "/";
CreateDirectory[targetDir, CreateIntermediateDirectories -> True];
Put[regionSummary, targetDir <> "Compare-RegionInfo-" <> famname <> ".m"];
Print["  Exported ", Length[regionSummary], " region(s) to verify/", famname, "/Compare-RegionInfo-", famname, ".m"];

(* Step 3: Export binary matrices *)
Print["\n[Step 3] Exporting binary matrices..."];
Alist = data["Family", "AList"];
ne = Length[Alist];
char = config["Modulus"];

binFile = targetDir <> "IBPMat_" <> famname <> ".bin";
ringFile = targetDir <> "RingData_" <> famname <> ".bin";

ExportIBPMatrixBinary`ExportBinaryIBPMatrix[binFile, expregdata, char];
ExportIBPMatrixBinary`ExportBinaryRingData[ringFile, expregdata, Alist, ne, char];

(* Export timing summary for verification log *)
mmaTiming = <|
    "FamilyDefinition" -> timeDefine,
    "RegionSolving" -> timeRegions
|>;
Export[targetDir <> "Compare-MMATiming-" <> famname <> ".m", mmaTiming, "Text"];

Print["  ", binFile, " + ", ringFile, " -> ", targetDir];

Print["\n========================================"];
Print["Compare-FamilyGenerate complete!"];
Print["  Next: ./build/test_expandFF ", famname];
Print["========================================"];
