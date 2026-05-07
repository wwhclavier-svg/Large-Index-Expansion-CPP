(* ::Package:: *)
(* Verify-Blade.wl — Blade IBP reduction verification of C++ LIE relations *)
(* Usage standalone: wolframscript -file Verify-Blade.wl <famname> [nuSample] *)
(* Load in-process:  $LoadedByBladeFamily=True; Get["Verify-Blade.wl"] *)

If[!TrueQ[$LoadedByBladeFamily],
    (* === STANDALONE MODE === *)
    current = If[$FrontEnd === Null, $InputFileName, NotebookFileName[]] // DirectoryName;
    $ProjectRoot = FileNameJoin[{current, ".."}];
    $BladePath = "/home/ykm/blade-workspace";

    Get[FileNameJoin[{$BladePath, "Blade.wl"}]];
    Get[FileNameJoin[{$ProjectRoot, "verify/FamilyDatabase/FamilyDatabase.wl"}]];

    If[Length[$ScriptCommandLine] >= 2, famname = $ScriptCommandLine[[2]], famname = "bub00"];
    nuSampleOverride = If[Length[$ScriptCommandLine] >= 3,
        ToExpression[$ScriptCommandLine[[3]]], Automatic];

    If[!KeyExistsQ[$FamilyDatabase, famname],
        Print["ERROR: Unknown family '", famname, "'"]; Exit[1]];

    config = $FamilyDatabase[famname];
    bladeName = StringReplace[famname, "-" -> ""];
    family = Symbol[bladeName];
    dimension = 4 - 2*eps;
    loop = config["LoopMomenta"]; leg = config["ExternalMomenta"];
    conservation = {}; replacement = config["KinematicRules"];
    propagator = config["Propagators"]; topsector = config["TopSector"];
    numeric = DeleteCases[config["Numeric"], "d" -> _];
    dValue = "d" /. config["Numeric"]; modulus = config["Modulus"];
    epsSol = Solve[dimension - dValue == 0, eps][[1]];

    Print["====================================================="];
    Print["  MODULE: BladeVerify  fam=", famname];
    Print["====================================================="];

    BLSetReducerOptions["MaxIncrement" -> 5, "BlackBoxRank" -> 3, "BlackBoxDot" -> 1];
    BLFamilyDefine[family, dimension, propagator, loop, leg,
        conservation, replacement, topsector, numeric];

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
    If[!ValueQ[$AllRelations], Print["ERROR: $AllRelations not defined"]; Exit[1]];

    NE = $AllRelations[[1]]["NE"];
    mkNuSample[] := If[nuSampleOverride =!= Automatic, nuSampleOverride,
        If[NE <= 2, {2, 5}, If[NE <= 5, {1, 2, 1, 2, 1}, Table[RandomInteger[{1, 15}], {NE}]]]];
];

(* === VERIFICATION LOGIC === *)
(* Shared cache — reuse if already present, create if not *)
If[!ValueQ[redCache], redCache = <||>];
reduceOne[e_List] := Module[{t, r}, t = BL[family, e]; r = BLReduce[{t}]; If[Length[r] > 0, r[[1]], 0]];
cachedReduce[e_List] := (If[!KeyExistsQ[redCache, e], redCache[e] = reduceOne[e]]; redCache[e]);

Print["\n=== MODULE 1: BladeVerify ============================="];

nuSample = mkNuSample[];
Print["  nu: ", nuSample, "  eps: ", eps /. epsSol];

Print["\n  [1/3] Top sector reduction..."];
topExps = ReplacePart[ConstantArray[0, Length[topsector]], Thread[Position[topsector, 1] -> 1]];
time1 = AbsoluteTiming[BLReduce[{BL[family, topExps]}];][[1]];
masters = BLFamilyInf[family]["Masters"];
Print["  Masters: ", Length[masters], "  (", Round[time1, 0.1], "s)"];
Do[Print["    ", m], {m, Take[masters, UpTo[8]]}];
If[Length[masters] > 8, Print["    ... and ", Length[masters]-8, " more"]];

Print["\n  [2/3] Reducing integrals..."];
allGexpSet = {};
Do[
    entry = $AllRelations[[i]];
    If[entry["NumSolutions"] == 0, Continue[]];
    alphas = entry["Alphas"];
    Do[AppendTo[allGexpSet, Table[nuSample[[j]] - alphas[[ai, j]], {j, NE}]], {ai, Length[alphas]}];
, {i, Length[$AllRelations]}];
allGexpSet = DeleteDuplicates[allGexpSet];
Print["  Unique integrals: ", Length[allGexpSet]];
time2 = AbsoluteTiming[Do[cachedReduce[g], {g, allGexpSet}];][[1]];
Print["  Reduced in ", Round[time2, 0.1], "s"];

Print["\n  [3/3] Verifying..."];
bladeTotalPass = 0; bladeTotalFail = 0;

Do[
    entry = $AllRelations[[i]];
    lev = entry["Lev"]; deg = entry["Deg"]; sols = entry["NumSolutions"];
    If[sols == 0, Continue[]];
    alphas = entry["Alphas"]; betas = entry["Betas"]; coeffMat = entry["Coefficients"];
    nAlpha = Length[alphas]; nBeta = Length[betas];

    passCount = 0; failCount = 0;
    Do[
        expr = 0; idx = 1;
        Do[
            alpha = alphas[[ai]];
            Do[
                beta = betas[[bi]]; coeff = coeffMat[[idx, si+1]];
                If[coeff != 0,
                    gExp = Table[nuSample[[j]] - alpha[[j]], {j, NE}];
                    betaFactor = Product[nuSample[[j]]^beta[[j]], {j, NE}];
                    expr = expr + coeff * betaFactor * cachedReduce[gExp];
                ];
                idx++;
            , {bi, nBeta}];
        , {ai, nAlpha}];

        expanded = Expand[expr /. epsSol];
        passQ = True;
        Do[
            coeff = Coefficient[expanded, m];
            If[coeff =!= 0 && NumericQ[coeff] && Mod[Numerator[coeff], modulus] =!= 0,
                passQ = False];
        , {m, masters}];
        If[passQ, passCount++, failCount++];
    , {si, 0, sols-1}];

    bladeTotalPass += passCount; bladeTotalFail += failCount;
    status = If[failCount == 0, "PASS", "FAIL (" <> ToString[failCount] <> "/" <> ToString[sols] <> ")"];
    Print["  (lev=", StringPadRight[ToString[lev], 3],
        " deg=", StringPadRight[ToString[deg], 3], ") ",
        StringPadRight[status, 18], " pass=", passCount, "/", sols];
, {i, Length[$AllRelations]}];

Print["\n  BladeVerify: ", bladeTotalPass, "/",
    bladeTotalPass + bladeTotalFail, " PASSED",
    If[bladeTotalFail > 0, " (" <> ToString[bladeTotalFail] <> " failed)", ""]];

If[!TrueQ[$LoadedByBladeFamily], Exit[If[bladeTotalFail > 0, 1, 0]]];
