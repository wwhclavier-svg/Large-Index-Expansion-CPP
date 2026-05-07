(* ::Package:: *)
(* Verify-Series.wl — Series expansion substitution verification *)
(* Usage standalone: wolframscript -file Verify-Series.wl <famname> *)
(* Load in-process:  $LoadedByBladeFamily=True; Get["Verify-Series.wl"] *)
(* Prerequisites: RelationMeta_<fam>.m + Compare-CPPResult-<fam>.m + ExpansionMMA_<fam>.m *)

If[!TrueQ[$LoadedByBladeFamily],
    current = If[$FrontEnd === Null, $InputFileName, NotebookFileName[]] // DirectoryName;
    $ProjectRoot = FileNameJoin[{current, ".."}];
    Get[FileNameJoin[{$ProjectRoot, "verify/FamilyDatabase/FamilyDatabase.wl"}]];

    If[Length[$ScriptCommandLine] >= 2, famname = $ScriptCommandLine[[2]], famname = "bub00"];
    If[!KeyExistsQ[$FamilyDatabase, famname],
        Print["ERROR: Unknown family '", famname, "'"]; Exit[1]];
    config = $FamilyDatabase[famname];
    modulus = config["Modulus"];

    allRelFiles = FileNames["AllRelations_" <> famname <> "_k*.m", $ProjectRoot];
    If[allRelFiles === {},
        Print["ERROR: No AllRelations_", famname, "_k*.m found"]; Exit[1]];
    allRelFilesNum = Select[allRelFiles,
        StringMatchQ[FileNameTake[#],
            "AllRelations_" <> famname <> "_k" ~~ DigitCharacter .. ~~ ".m"] &];
    extractK[fname_String] := ToExpression[
        StringCases[fname, "AllRelations_" <> famname <> "_k" ~~ k : DigitCharacter .. ~~ ".m" :> k][[1]]];
    allRelFile = Last[SortBy[allRelFilesNum, extractK]];
    Get[allRelFile];
    NE = $AllRelations[[1]]["NE"];
];

Print["\n=== MODULE 2: SeriesVerify ============================"];

(* ====================== Load expansion coefficients ====================== *)
$ExpansionFile = FileNameJoin[{$ProjectRoot, "ExpansionMMA_" <> famname <> ".m"}];
$CPPResultFile = FileNameJoin[{$ProjectRoot, "verify", famname, "Compare-CPPResult-" <> famname <> ".m"}];

cResultFile = If[FileExistsQ[$CPPResultFile], $CPPResultFile,
    If[FileExistsQ[$ExpansionFile], $ExpansionFile,
        Print["  SKIP: No expansion file found (ExpansionMMA or Compare-CPPResult)"];
        If[!TrueQ[$LoadedByBladeFamily], Exit[0], Abort[]]
    ]];
Get[cResultFile];
If[!ValueQ[$ExpansionResults],
    Print["  SKIP: $ExpansionResults not defined in ", cResultFile];
    If[!TrueQ[$LoadedByBladeFamily], Exit[0], Abort[]]];

(* ====================== Load A_i values from RelationMeta ====================== *)
$MetaFile = FileNameJoin[{$ProjectRoot, "RelationMeta_" <> famname <> ".m"}];
If[FileExistsQ[$MetaFile],
    Get[$MetaFile];
    If[ValueQ[$RelationMeta],
        regimes = $RelationMeta["Regimes"];
        topRegime = regimes[[1]];
        aVals = topRegime["A"];
        aInvVals = topRegime["Ainv"];
        Print["  A-values from RelationMeta: ", aVals];
        ,
        Print["  WARNING: $RelationMeta not defined, using hardcoded A=59808223"];
        aVals = Table[59808223, {NE}];
        aInvVals = aVals;
    ];
,
    Print["  WARNING: RelationMeta not found, using hardcoded A=59808223"];
    aVals = Table[59808223, {NE}];
    aInvVals = aVals;
];

(* ====================== Verification ====================== *)
seriesTotalPass = 0; seriesTotalFail = 0;
order = 4;
secpos = topsector /. {0 -> 0, 1 -> 1};
vlistSym = Table[Symbol["v" <> ToString[i]], {i, NE}];
varRule = Table[A[i] -> aInvVals[[i]], {i, NE}];

Do[
    entry = $AllRelations[[i]];
    lev = entry["Lev"]; deg = entry["Deg"]; sols = entry["NumSolutions"];
    If[sols == 0, Continue[]];
    alphas = entry["Alphas"]; betas = entry["Betas"];
    coeffMat = entry["Coefficients"];

    hlist = None;
    Do[
        epEntry = $ExpansionResults[[j, 1]];
        If[KeyExistsQ[epEntry, "Lev"] && KeyExistsQ[epEntry, "Deg"],
            If[epEntry["Lev"] == lev && epEntry["Deg"] == deg,
                hlist = epEntry["Solutions"][[1]]["H"]; Break[]],
            If[hlist === None, hlist = epEntry["Solutions"][[1]]["H"]]
        ];
    , {j, Length[$ExpansionResults]}];
    If[hlist === None,
        Print["  (lev=", lev, ", deg=", deg, ") SKIP: no expansion data"];
        Continue[]];

    hExpr = Sum[1/n^k * (hlist[[k+1]] /. Thread[vlistSym -> vlistSym]),
                  {k, 0, Length[hlist]-1}];

    relExpr = 0; idx = 1;
    Do[
        alpha = alphas[[ai]];
        Do[
            beta = betas[[bi]];
            coeff = coeffMat[[idx, 1]];
            If[coeff != 0,
                betaPoly = Times @@ Thread[Power[vlistSym, beta]];
                gShift = vlistSym - Table[alpha[[j]], {j, NE}];
                shiftExpr = PolynomialMod[
                    Times @@ Thread[Power[Alist /. varRule, -gShift]] *
                    (hExpr /. Thread[vlistSym -> gShift]), modulus];
                relExpr = relExpr + coeff * betaPoly * shiftExpr;
            ];
            idx++;
        , {bi, Length[betas]}];
    , {ai, Length[alphas]}];

    relSeries = (relExpr / n^deg) /. Thread[vlistSym -> vlistSym + secpos*n];
    relSeries = CoefficientList[relSeries, 1/n, order+1];
    relSeries = PolynomialMod[#, modulus] & /@ Flatten[{relSeries}];

    passQ = AllTrue[relSeries, # === 0 &];
    If[passQ, seriesTotalPass++, seriesTotalFail++];
    Print["  (lev=", lev, ", deg=", deg,
        ") Series: ", If[passQ, "PASS", "FAIL"]];
, {i, Length[$AllRelations]}];

Print["\n  SeriesVerify: ",
    If[seriesTotalPass + seriesTotalFail > 0,
        ToString[seriesTotalPass] <> "/" <> ToString[seriesTotalPass + seriesTotalFail] <>
        If[seriesTotalFail > 0,
            " (" <> ToString[seriesTotalFail] <> " failed)", " PASS"],
        "N/A"]];

If[!TrueQ[$LoadedByBladeFamily], Exit[If[seriesTotalFail > 0, 1, 0]]];
