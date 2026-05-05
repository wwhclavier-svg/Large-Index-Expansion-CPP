(* ::Package:: *)
(* Compare-VerifyLog.wl — Unified Verification Log Generator *)
(* Usage: wolframscript -file Compare-VerifyLog.wl <famname> *)
(* Reads all intermediate files and generates VerifyLog-<famname>-<date>.md *)

SetDirectory[DirectoryName[$InputFileName]];

(* Parse CLI argument *)
If[Length[$ScriptCommandLine] < 2,
    Print["ERROR: Missing family name argument."];
    Print["Usage: wolframscript -file Compare-VerifyLog.wl <famname>"];
    Exit[1]
];
famname = $ScriptCommandLine[[2]];

modulus = 179424673; (* Prime[10000000] *)
dateStr = DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}];

Print["========================================"];
Print["Compare-VerifyLog     fam = ", famname];
Print["  Date = ", dateStr];
Print["========================================"];

(* Search for a file across multiple base directories *)
findFile[baseName_] := Module[{searchPaths, found},
    searchPaths = {"",
                   "../../",
                   "../../build/",
                   "../" <> famname <> "/"};
    found = SelectFirst[searchPaths, FileExistsQ[# <> baseName] &];
    If[MissingQ[found], baseName, found <> baseName]
];

(* ---- Load files ---- *)

(* 1. Region info (from Compare-FamilyGenerate.wl) *)
regionFile = findFile["Compare-RegionInfo-" <> famname <> ".m"];
regionSummary = If[FileExistsQ[regionFile],
    Get[regionFile],
    Print["[WARN] Region info not found: ", regionFile];
    {}
];

(* 2. C++ metadata *)
cppMetaFile = findFile["Compare-CPPMeta-" <> famname <> ".m"];
cppMeta = If[FileExistsQ[cppMetaFile],
    Get[cppMetaFile];
    If[ValueQ[$CPPMeta], $CPPMeta, <||>],
    Print["[WARN] C++ metadata not found: ", cppMetaFile];
    <||>
];

(* 3. MMA timing from Compare-FamilyGenerate *)
mmaTimingFile = findFile["Compare-MMATiming-" <> famname <> ".m"];
mmaTimingData = If[FileExistsQ[mmaTimingFile],
    Get[mmaTimingFile],
    Print["[WARN] MMA timing not found: ", mmaTimingFile];
    <||>
];

(* 4. MMA expansion timing from Compare-Expand *)
mmaExpandFile = findFile["Compare-MMAExpandTiming-" <> famname <> ".m"];
mmaExpandData = If[FileExistsQ[mmaExpandFile],
    Get[mmaExpandFile],
    Print["[WARN] MMA expansion timing not found: ", mmaExpandFile];
    <||>
];

(* 5. C++ expansion results *)
cppResultFile = findFile["Compare-CPPResult-" <> famname <> ".m"];
cHlist = If[FileExistsQ[cppResultFile],
    Get[cppResultFile];
    $ExpansionResults[[1, 1, "Solutions", 1, "H"]],
    Print["[WARN] C++ results not found: ", cppResultFile];
    {}
];

(* 6. MMA expansion results *)
mmaResultFile = findFile["Compare-MMAResult-" <> famname <> ".m"];
mHlist = If[FileExistsQ[mmaResultFile],
    Get[mmaResultFile];
    $ExpansionResults[[1, 1, "Solutions", 1, "H"]],
    Print["[WARN] MMA results not found: ", mmaResultFile];
    {}
];

(* ---- Helper: format solution text ---- *)
formatSolution[rule_] := If[Length[rule] === 0,
    "(none -- no parametrized variables)",
    StringRiffle[rule /. Rule[a_, b_] :> ToString[a] <> " -> " <> ToString[b, InputForm], ", "]
];

(* ---- Helper: compress polynomial to shorter form ---- *)
compressPoly[poly_, maxTerms_:5] := Module[{terms, nTerms},
    If[poly === 0, Return["0"]];
    terms = If[Head[poly] === Plus, List @@ poly, {poly}];
    nTerms = Length[terms];
    If[nTerms <= maxTerms,
        ToString[poly, InputForm],
        ToString[Plus @@ Take[terms, maxTerms], InputForm] <> " + ... (" <> ToString[nTerms] <> " terms)"
    ]
];

(* ---- Build Markdown ---- *)
targetDir = "../" <> famname <> "/";
CreateDirectory[targetDir, CreateIntermediateDirectories -> True];
md = OpenWrite[targetDir <> "VerifyLog-" <> famname <> ".md", PageWidth -> Infinity];

WriteLine[md, "# Verify Log: " <> famname];
WriteLine[md, "**Date:** " <> dateStr];
WriteLine[md, "**Modulus:** " <> ToString[modulus]];

(* Try to get description from FamilyDatabase *)
desc = "N/A";
famdbFile = "../FamilyDatabase/FamilyDatabase.wl";
If[FileExistsQ[famdbFile],
    Get[famdbFile];
    If[KeyExistsQ[$FamilyDatabase, famname],
        desc = $FamilyDatabase[famname]["Description"]
    ]
];
WriteLine[md, "**Description:** " <> desc];
WriteLine[md, ""];
WriteLine[md, "---"];
WriteLine[md, ""];

(* ================================================================ *)
(* Part 1: Characteristic Equation *)
(* ================================================================ *)
WriteLine[md, "## 1. Characteristic Equation: Minimal Associated Primes"];
WriteLine[md, ""];

If[Length[regionSummary] === 0,
    WriteLine[md, "*(No region data available -- run Compare-FamilyGenerate.wl first)*"];
    WriteLine[md, ""];
,
    (* Group by sector *)
    sectors = Union[#Sector & /@ regionSummary];

    (* Sector summary table *)
    WriteLine[md, "### Summary"];
    WriteLine[md, ""];
    WriteLine[md, "| Sector | #Regions | Total #Solutions |"];
    WriteLine[md, "|--------|:--------:|:----------------:|"];
    Do[
        regsInSector = Select[regionSummary, #Sector === sector &];
        totalSol = Total[#NumSolutions & /@ regsInSector];
        nRegs = Length[regsInSector];
        WriteLine[md, "| " <> ToString[sector, InputForm] <> " | " <> ToString[nRegs] <> " | " <> ToString[totalSol] <> " |"];
    , {sector, sectors}];
    WriteLine[md, ""];

    (* Detail per sector *)
    Do[
        regsInSector = Select[regionSummary, #Sector === sector &];
        WriteLine[md, "### Sector " <> ToString[sector, InputForm]];
        WriteLine[md, ""];

        Do[
            reg = regsInSector[[ri]];
            WriteLine[md, "#### Primary Component (" <> ToString[ri] <> ")"];

            (* Leading degree from VarDeg *)
            leadDeg = If[Length[reg["VarDeg"]] > 0,
                StringRiffle[
                    Table["`A[" <> ToString[i] <> "], " <> ToString[reg["VarDeg"][[i]]] <> "`",
                          {i, Length[reg["VarDeg"]]}],
                    ", "
                ],
                "`{}`"
            ];

            numSol = Times @@ reg["VarDeg"];
            varDegStr = StringRiffle[ToString /@ reg["VarDeg"], "x"];

            WriteLine[md, ""];
            WriteLine[md, "| Field | Value |"];
            WriteLine[md, "|-------|-------|"];
            WriteLine[md, "| Leading Degree | `{{" <> leadDeg <> "}}` |"];
            WriteLine[md, "| Generator (indep) | `" <> ToString[reg["VarIndep"], InputForm] <> "` |"];
            WriteLine[md, "| Parametrized (dep) | `" <> ToString[reg["VarDep"], InputForm] <> "` |"];
            WriteLine[md, "| **#Solutions** | **" <> ToString[numSol] <> "** (= " <> varDegStr <> ") |"];
            WriteLine[md, "| Solution | `" <> formatSolution[reg["VarRule"]] <> "` |"];
            WriteLine[md, "| Min.Poly. | `" <> If[Length[reg["MinPoly"]] === 0 || (ListQ[reg["MinPoly"]] && Length[reg["MinPoly"][[1]]] === 0),
                "(none -- rational)", ToString[reg["MinPoly"], InputForm]] <> "` |"];
            WriteLine[md, ""];
        , {ri, Length[regsInSector]}];
    , {sector, sectors}];
];

(* ================================================================ *)
(* Part 2: Determined Parameters *)
(* ================================================================ *)
WriteLine[md, "## 2. Determined Parameters"];
WriteLine[md, ""];

cRegions = Lookup[cppMeta, "Regions", {}];
firstRegion = If[Length[cRegions] > 0, First[cRegions], <||>];

increVal = Lookup[firstRegion, "Incre", "?"];
nimaxVal = Lookup[firstRegion, "Nimax", "?"];
neVal = Lookup[firstRegion, "NE", "?"];
nbVal = Lookup[firstRegion, "NB", "?"];

WriteLine[md, "| Parameter | Value | Description |"];
WriteLine[md, "|-----------|-------|-------------|"];
WriteLine[md, "| `incre` | " <> ToString[increVal] <> " | nu-polynomial expansion growth per order |"];
WriteLine[md, "| `nimax` | " <> ToString[nimaxVal] <> " | Solution space dimension (independent 1/n series) |"];
WriteLine[md, "| `ne` | " <> ToString[neVal] <> " | #external variables (A_i) |"];
WriteLine[md, "| `nb` | " <> ToString[nbVal] <> " | Monomial basis dimension (coordinate ring) |"];

(* Also show per-region breakdown if multiple regions *)
If[Length[cRegions] > 1,
    WriteLine[md, ""];
    WriteLine[md, "### Per-Region Details"];
    WriteLine[md, ""];
    WriteLine[md, "| Region | NE | NB | NIBP | Incre | Nimax |"];
    WriteLine[md, "|--------|:--:|:--:|:----:|:-----:|:-----:|"];
    Do[
        r = cRegions[[i]];
        WriteLine[md, "| " <> ToString[Lookup[r, "RegionIndex", i]] <>
                  " | " <> ToString[Lookup[r, "NE", "?"]] <>
                  " | " <> ToString[Lookup[r, "NB", "?"]] <>
                  " | " <> ToString[Lookup[r, "NIBP", "?"]] <>
                  " | " <> ToString[Lookup[r, "Incre", "?"]] <>
                  " | " <> ToString[Lookup[r, "Nimax", "?"]] <> " |"];
    , {i, Length[cRegions]}];
];
WriteLine[md, ""];

(* ================================================================ *)
(* Part 3: Computation Time *)
(* ================================================================ *)
WriteLine[md, "## 3. Computation Time"];
WriteLine[md, ""];

mmaDefineTime = Lookup[mmaTimingData, "FamilyDefinition", "?"];
mmaRegionTime = Lookup[mmaTimingData, "RegionSolving", "?"];
mmaExpandTime = Lookup[mmaExpandData, "Expansion", "?"];
cppTotalTime = Lookup[cppMeta, "TotalTime", "?"];

mmaTotal = If[NumericQ[mmaDefineTime] && NumericQ[mmaRegionTime] && NumericQ[mmaExpandTime],
    mmaDefineTime + mmaRegionTime + mmaExpandTime,
    "?"
];

WriteLine[md, "| Step | MMA (s) | C++ (s) |"];
WriteLine[md, "|------|:------:|:------:|"];
WriteLine[md, "| Family Definition | " <> If[NumericQ[mmaDefineTime], ToString[NumberForm[mmaDefineTime, {4, 3}]], ToString[mmaDefineTime]] <> " | -- |"];
WriteLine[md, "| Region Solving (GB + minAssGTZ) | " <> If[NumericQ[mmaRegionTime], ToString[NumberForm[mmaRegionTime, {4, 3}]], ToString[mmaRegionTime]] <> " | -- |"];
WriteLine[md, "| Expansion | " <> If[NumericQ[mmaExpandTime], ToString[NumberForm[mmaExpandTime, {4, 3}]], ToString[mmaExpandTime]] <> " | " <> If[NumericQ[cppTotalTime], ToString[NumberForm[cppTotalTime, {6, 6}]], ToString[cppTotalTime]] <> " |"];
WriteLine[md, "| **Total** | **" <> If[NumericQ[mmaTotal], ToString[NumberForm[mmaTotal, {4, 3}]], ToString[mmaTotal]] <> "** | **" <> If[NumericQ[cppTotalTime], ToString[NumberForm[cppTotalTime, {6, 6}]], ToString[cppTotalTime]] <> "** |"];
WriteLine[md, ""];

(* ================================================================ *)
(* Part 4: Series Expansion — Orders 0 & 1 *)
(* ================================================================ *)
WriteLine[md, "## 4. Series Expansion: Orders 0 and 1"];
WriteLine[md, ""];

maxDetail = Min[1, Min[If[Length[cHlist] > 0, Length[cHlist] - 1, -1],
                       If[Length[mHlist] > 0, Length[mHlist] - 1, -1]]];

For[k = 0, k <= maxDetail, k++,
    WriteLine[md, "### Order k=" <> ToString[k]];
    WriteLine[md, ""];

    cppPoly = If[Length[cHlist] > k, cHlist[[k + 1]], Missing[]];
    mmaPoly = If[Length[mHlist] > k, mHlist[[k + 1]], Missing[]];

    If[MissingQ[cppPoly] || MissingQ[mmaPoly],
        WriteLine[md, "*(data not available)*"];
        WriteLine[md, ""];
        Continue[]
    ];

    diff = PolynomialMod[Expand[mmaPoly - cppPoly], modulus];
    status = If[diff === 0, "[MATCH]", "[MISMATCH]"];

    WriteLine[md, "| | Expression |"];
    WriteLine[md, "|---|-----------|"];
    WriteLine[md, "| MMA | `" <> ToString[mmaPoly, InputForm] <> "` |"];
    WriteLine[md, "| C++ | `" <> ToString[cppPoly, InputForm] <> "` |"];
    WriteLine[md, "| Diff | `" <> ToString[diff, InputForm] <> "` **" <> status <> "** |"];
    WriteLine[md, ""];
];

(* ================================================================ *)
(* Part 5: Full Comparison Summary *)
(* ================================================================ *)
WriteLine[md, "## 5. Full Comparison Summary"];
WriteLine[md, ""];

maxFull = Min[If[Length[cHlist] > 0, Length[cHlist] - 1, -1],
              If[Length[mHlist] > 0, Length[mHlist] - 1, -1]];

WriteLine[md, "| k | MMA | C++ | Diff | Result |"];
WriteLine[md, "|:--|------|------|------|:------:|"];

allMatch = True;
For[k = 0, k <= maxFull, k++,
    cppPoly = cHlist[[k + 1]];
    mmaPoly = mHlist[[k + 1]];
    diff = PolynomialMod[Expand[mmaPoly - cppPoly], modulus];
    result = If[diff === 0, "MATCH", "MISMATCH"];
    If[diff =!= 0, allMatch = False];

    (* Shorten polynomials for table *)
    cppShort = compressPoly[cppPoly, 3];
    mmaShort = compressPoly[mmaPoly, 3];
    diffShort = compressPoly[diff, 3];

    WriteLine[md, "| " <> ToString[k] <> " | `" <> cppShort <> "` | `" <> mmaShort <> "` | `" <> diffShort <> "` | " <> result <> " |"];
];
WriteLine[md, ""];

verdict = If[allMatch, "[PASS] -- All " <> ToString[maxFull + 1] <> " orders match.", "[FAIL] -- Some orders do not match. See details above."];
WriteLine[md, "**Verdict: " <> verdict <> "**"];
WriteLine[md, ""];

(* ---- Write file ---- *)
Close[md];

logFile = targetDir <> "VerifyLog-" <> famname <> ".md";
Print[""];
Print["========================================"];
Print["Verification log generated: ", logFile];
If[allMatch,
    Print["Verdict: [PASS] All orders match!"],
    Print["Verdict: [FAIL] Some mismatches detected."]
];
Print["========================================"];

If[!allMatch, Exit[1]];
