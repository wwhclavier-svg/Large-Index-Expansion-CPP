(* ::Package:: *)
(* VerifyExpand-Compare.wl — C++ vs MMA comparison + verification log *)
(* Usage: wolframscript -file VerifyExpand-Compare.wl <famname> *)
(*
    Reads:
      C++  result:  verify/<fam>/Compare-CPPResult-<fam>.m
      MMA  result:  verify/<fam>/VerifyExpansion-MMAExpansion.m
    Generates:
      verify/<fam>/VerifyLog-<fam>.md
*)

SetDirectory[DirectoryName[$InputFileName]];
Get["VerifyExpand-Prepare.wl"];  (* loads shared utilities *)

(* ---- Parse CLI ---- *)
famname = parseFamilyArg[];

modulus = 179424673;  (* Prime[10000000] *)
dateStr = DateString[{"Year", "-", "Month", "-", "Day", " ", "Hour", ":", "Minute", ":", "Second"}];
target = targetDir[famname];

Print["========================================"];
Print["VerifyExpand-Compare     fam = ", famname];
Print["  Modulus = ", modulus];
Print["  Date = ", dateStr];
Print["========================================"];

(* ============================================================ *)
(* SECTION: Load results                                        *)
(* ============================================================ *)

(* 1. C++ result *)
cppBase = "Compare-CPPResult-" <> famname <> ".m";
cppFile = findResultFile[famname, cppBase];
Print["\n[Step 1] Loading C++ result from: ", cppFile];
If[!FileExistsQ[cppFile],
    Print["[FAIL] C++ result not found."];
    Print["  Run: cd ../.. && ./build/test_expandFF ", famname];
    Exit[1]
];
Get[cppFile];
hlistCPP = $ExpansionResults[[1, 1, "Solutions", 1, "H"]];
Print["  C++ result loaded. Orders: ", Length[hlistCPP] - 1];

(* 2. C++ metadata (for log) *)
cppMetaFile = findResultFile[famname, "Compare-CPPMeta-" <> famname <> ".m"];
cppMeta = If[FileExistsQ[cppMetaFile],
    Get[cppMetaFile];
    If[ValueQ[$CPPMeta], $CPPMeta, <||>],
    Print["[WARN] C++ metadata not found"];
    <||>
];

(* 3. MMA result *)
mmaBase = "VerifyExpansion-MMAExpansion.m";
mmaFile = findResultFile[famname, mmaBase];
Print["\n[Step 2] Loading MMA result from: ", mmaFile];
If[!FileExistsQ[mmaFile],
    Print["[FAIL] MMA result not found."];
    Print["  Run: wolframscript -file VerifyExpand-MMAExpand.wl ", famname];
    Exit[1]
];
Get[mmaFile];
hlistMMA = $ExpansionResults[[1, 1, "Solutions", 1, "H"]];
Print["  MMA result loaded. Orders: ", Length[hlistMMA] - 1];

(* 4. Region info *)
regionFile = findResultFile[famname, "Compare-RegionInfo-" <> famname <> ".m"];
regionSummary = If[FileExistsQ[regionFile],
    Get[regionFile],
    Print["[WARN] Region info not found"];
    {}
];

(* 5. Timing data *)
mmaTimingFile = findResultFile[famname, "Compare-MMATiming-" <> famname <> ".m"];
mmaTimingData = If[FileExistsQ[mmaTimingFile], Get[mmaTimingFile], <||>];

mmaExpandFile = findResultFile[famname, "Compare-MMAExpandTiming-" <> famname <> ".m"];
mmaExpandData = If[FileExistsQ[mmaExpandFile], Get[mmaExpandFile], <||>];

(* ============================================================ *)
(* SECTION: Compare coefficients                                *)
(* ============================================================ *)
Print["\n[Step 3] Comparing coefficients..."];

maxk = Min[Min[Length[hlistCPP], Length[hlistMMA]] - 1, 4];
allMatch = True;

For[k = 0, k <= maxk, k++,
    cppPoly = hlistCPP[[k + 1]];
    mmaPoly = hlistMMA[[k + 1]];
    diff = PolynomialMod[Expand[mmaPoly - cppPoly], modulus];

    If[diff === 0,
        Print["  [MATCH] k=", k],
        Print["  [MISMATCH] k=", k];
        Print["    C++ = ", cppPoly];
        Print["    MMA = ", mmaPoly];
        Print["    diff = ", diff];
        allMatch = False;
    ];
];

If[maxk < 4,
    Print["  (only orders 0..", maxk, " available)"]
];

(* ============================================================ *)
(* SECTION: Generate verification log                           *)
(* ============================================================ *)
Print["\n[Step 4] Generating verification log..."];

md = OpenWrite[target <> "VerifyLog-" <> famname <> ".md", PageWidth -> Infinity];

moduleLogHeader[] := (
    WriteLine[md, "# Verify Log: " <> famname];
    WriteLine[md, "**Date:** " <> dateStr];
    WriteLine[md, "**Modulus:** " <> ToString[modulus]];

    Quiet[Get[$CPPVerifyRoot <> "FamilyDatabase/FamilyDatabase.wl"]];
    If[ValueQ[$FamilyDatabase] && KeyExistsQ[$FamilyDatabase, famname],
        WriteLine[md, "**Description:** " <> $FamilyDatabase[famname]["Description"]];
    ];
    WriteLine[md, ""];
    WriteLine[md, "---"];
    WriteLine[md, ""];
);

moduleRegionSummary[] := (
    WriteLine[md, "## 1. Characteristic Equation: Minimal Associated Primes"];
    WriteLine[md, ""];
    If[Length[regionSummary] === 0,
        WriteLine[md, "*(No region data — run VerifyExpand-Prepare.wl first)*"];
        WriteLine[md, ""];
        Return[];
    ];

    sectors = Union[#Sector & /@ regionSummary];

    (* Sector summary table *)
    WriteLine[md, "### Summary"];
    WriteLine[md, ""];
    WriteLine[md, "| Sector | #Regions | Total #Solutions |"];
    WriteLine[md, "|--------|:--------:|:----------------:|"];
    Do[
        regsInSec = Select[regionSummary, #Sector === sector &];
        totalSol = Total[#NumSolutions & /@ regsInSec];
        WriteLine[md, "| " <> ToString[sector, InputForm] <> " | " <>
            ToString[Length[regsInSec]] <> " | " <> ToString[totalSol] <> " |"];
    , {sector, sectors}];
    WriteLine[md, ""];

    (* Detail per sector *)
    Do[
        regsInSec = Select[regionSummary, #Sector === sector &];
        WriteLine[md, "### Sector " <> ToString[sector, InputForm]];
        WriteLine[md, ""];
        Do[
            reg = regsInSec[[ri]];
            leadDeg = If[Length[reg["VarDeg"]] > 0,
                StringRiffle[
                    Table["`A[" <> ToString[i] <> "], " <> ToString[reg["VarDeg"][[i]]] <> "`", {i, Length[reg["VarDeg"]]}],
                    ", "
                ],
                "`{}`"
            ];
            WriteLine[md, ""];
            WriteLine[md, "| Field | Value |"];
            WriteLine[md, "|-------|-------|"];
            WriteLine[md, "| Leading Degree | `{{" <> leadDeg <> "}}` |"];
            WriteLine[md, "| Generator (indep) | `" <> ToString[reg["VarIndep"], InputForm] <> "` |"];
            WriteLine[md, "| Parametrized (dep) | `" <> ToString[reg["VarDep"], InputForm] <> "` |"];
            WriteLine[md, "| **#Solutions** | **" <> ToString[Times @@ reg["VarDeg"]] <> "** |"];
            WriteLine[md, "| Solution | `" <> formatSolution[reg["VarRule"]] <> "` |"];

            minPolyStr = If[Length[reg["MinPoly"]] === 0 || (ListQ[reg["MinPoly"]] && Length[reg["MinPoly"][[1]]] === 0),
                "(none -- rational)", ToString[reg["MinPoly"], InputForm]];
            WriteLine[md, "| Min.Poly. | `" <> minPolyStr <> "` |"];
            WriteLine[md, ""];
        , {ri, Length[regsInSec]}];
    , {sector, sectors}];
);

moduleParameters[] := (
    WriteLine[md, "## 2. Determined Parameters"];
    WriteLine[md, ""];
    cRegions = Lookup[cppMeta, "Regions", {}];
    firstRegion = If[Length[cRegions] > 0, First[cRegions], <||>];
    WriteLine[md, "| Parameter | Value | Description |"];
    WriteLine[md, "|-----------|-------|-------------|"];
    WriteLine[md, "| `incre` | " <> ToString[Lookup[firstRegion, "Incre", "?"]] <> " | growth per order |"];
    WriteLine[md, "| `nimax` | " <> ToString[Lookup[firstRegion, "Nimax", "?"]] <> " | solution space dim |"];
    WriteLine[md, "| `ne` | " <> ToString[Lookup[firstRegion, "NE", "?"]] <> " | #external variables |"];
    WriteLine[md, "| `nb` | " <> ToString[Lookup[firstRegion, "NB", "?"]] <> " | monomial basis dim |"];
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
);

moduleTiming[] := (
    WriteLine[md, "## 3. Computation Time"];
    WriteLine[md, ""];
    mmd = Lookup[mmaTimingData, "FamilyDefinition", "?"];
    mmr = Lookup[mmaTimingData, "RegionSolving", "?"];
    mme = Lookup[mmaExpandData, "Expansion", "?"];
    cpt = Lookup[cppMeta, "TotalTime", "?"];
    mmaTotal = If[NumericQ[mmd] && NumericQ[mmr] && NumericQ[mme], mmd + mmr + mme, "?"];
    WriteLine[md, "| Step | MMA (s) | C++ (s) |"];
    WriteLine[md, "|------|:------:|:------:|"];
    WriteLine[md, "| Family Definition | " <> If[NumericQ[mmd], ToString[NumberForm[mmd, {4,3}]], ToString[mmd]] <> " | -- |"];
    WriteLine[md, "| Region Solving | " <> If[NumericQ[mmr], ToString[NumberForm[mmr, {4,3}]], ToString[mmr]] <> " | -- |"];
    WriteLine[md, "| Expansion | " <> If[NumericQ[mme], ToString[NumberForm[mme, {4,3}]], ToString[mme]] <> " | " <> If[NumericQ[cpt], ToString[NumberForm[cpt, {6,6}]], ToString[cpt]] <> " |"];
    WriteLine[md, "| **Total** | **" <> If[NumericQ[mmaTotal], ToString[NumberForm[mmaTotal, {4,3}]], ToString[mmaTotal]] <> "** | **" <> If[NumericQ[cpt], ToString[NumberForm[cpt, {6,6}]], ToString[cpt]] <> "** |"];
    WriteLine[md, ""];
);

moduleDetailOrders[] := (
    WriteLine[md, "## 4. Series Expansion: Orders 0 and 1"];
    WriteLine[md, ""];
    maxDetail = Min[1, Min[Length[hlistCPP], Length[hlistMMA]] - 1];
    For[k = 0, k <= maxDetail, k++,
        WriteLine[md, "### Order k=" <> ToString[k]];
        WriteLine[md, ""];
        cppP = hlistCPP[[k + 1]]; mmaP = hlistMMA[[k + 1]];
        d = PolynomialMod[Expand[mmaP - cppP], modulus];
        status = If[d === 0, "[MATCH]", "[MISMATCH]"];
        WriteLine[md, "| | Expression |"];
        WriteLine[md, "|---|-----------|"];
        WriteLine[md, "| MMA | `" <> ToString[mmaP, InputForm] <> "` |"];
        WriteLine[md, "| C++ | `" <> ToString[cppP, InputForm] <> "` |"];
        WriteLine[md, "| Diff | `" <> ToString[d, InputForm] <> "` **" <> status <> "** |"];
        WriteLine[md, ""];
    ];
);

moduleSummaryTable[] := (
    WriteLine[md, "## 5. Full Comparison Summary"];
    WriteLine[md, ""];
    maxFull = Min[Length[hlistCPP], Length[hlistMMA]] - 1;
    WriteLine[md, "| k | MMA | C++ | Diff | Result |"];
    WriteLine[md, "|:--|------|------|------|:------:|"];
    Do[
        cppP = hlistCPP[[k + 1]]; mmaP = hlistMMA[[k + 1]];
        d = PolynomialMod[Expand[mmaP - cppP], modulus];
        result = If[d === 0, "MATCH", "MISMATCH"];
        If[d =!= 0, allMatch = False];
        WriteLine[md, "| " <> ToString[k] <> " | `" <> compressPoly[cppP, 3] <> "` | `" <> compressPoly[mmaP, 3] <> "` | `" <> compressPoly[d, 3] <> "` | " <> result <> " |"];
    , {k, 0, maxFull}];
    WriteLine[md, ""];
    verdict = If[allMatch, "[PASS] -- All " <> ToString[maxFull + 1] <> " orders match.", "[FAIL] -- Some orders do not match. See details above."];
    WriteLine[md, "**Verdict: " <> verdict <> "**"];
    WriteLine[md, ""];
);

(* ---- Build log ---- *)
moduleLogHeader[];
moduleRegionSummary[];
moduleParameters[];
moduleTiming[];
moduleDetailOrders[];
moduleSummaryTable[];

Close[md];

logFile = target <> "VerifyLog-" <> famname <> ".md";
Print["  Verdict: ", If[allMatch, "[PASS]", "[FAIL]"]];
Print["  Log: ", logFile];

Print["\n========================================"];
If[allMatch,
    Print["[PASS] All orders match! C++ and MMA are consistent."],
    Print["[FAIL] Some orders do not match. See details above."]
];
Print["========================================"];

If[!allMatch, Exit[1]];
