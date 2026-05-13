(* ::Package:: *)
(* Verify-Kira.wl — Kira IBP rule cross-validation of C++ LIE relations *)
(* Usage standalone: wolframscript -file Verify-Kira.wl <famname> *)
(* Load in-process:  $LoadedByBladeFamily=True; Get["Verify-Kira.wl"] *)

If[!TrueQ[$LoadedByBladeFamily],
    current = If[$FrontEnd === Null, $InputFileName, NotebookFileName[]] // DirectoryName;
    $ProjectRoot = FileNameJoin[{current, ".."}];
    Get[FileNameJoin[{$ProjectRoot, "families/FamilyDatabase.wl"}]];

    If[Length[$ScriptCommandLine] >= 2, famname = $ScriptCommandLine[[2]], famname = "bub00"];
    If[!KeyExistsQ[$FamilyDatabase, famname],
        Print["ERROR: Unknown family '", famname, "'"]; Exit[1]];

    config = $FamilyDatabase[famname];
    modulus = config["Modulus"];
    dValue = "d" /. config["Numeric"];

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

    mkNuSample[] := If[NE <= 2, {2, 5}, If[NE <= 5, {1, 2, 1, 2, 1}, Table[RandomInteger[{1, 15}], {NE}]]];
];

Print["\n=== MODULE 3: KiraVerify =============================="];

$KiraRulesFile = FileNameJoin[{$ProjectRoot, "verify", famname, "kira_integrals.m"}];
$KiraLoaderFile = FileNameJoin[{$ProjectRoot, "verify", "VerifyUtility", "KiraRuleLoader.wl"}];

kiraTotalPass = 0; kiraTotalFail = 0;

If[!FileExistsQ[$KiraRulesFile] || !FileExistsQ[$KiraLoaderFile],
    Print["  SKIP: Kira rules or KiraRuleLoader.wl not found"];
,
    Get[$KiraLoaderFile];
    If[!ValueQ[LoadKiraRules],
        Print["  SKIP: LoadKiraRules not defined in KiraRuleLoader.wl"];
    ,
        {allRulesSymbolic, kiraReduce} = LoadKiraRules[
            ParentDirectory[$KiraRulesFile],
            "d" -> dValue, "Modulus" -> modulus, "NProp" -> NE];
        Print["  Rules loaded: ", Length[allRulesSymbolic], " symbolic rules"];

        nuSampleKira = mkNuSample[];
        Do[
            entry = $AllRelations[[i]];
            lev = entry["Lev"]; deg = entry["Deg"]; sols = entry["NumSolutions"];
            If[sols == 0, Continue[]];
            alphas = entry["Alphas"]; betas = entry["Betas"];
            coeffMat = entry["Coefficients"];

            passCount = 0; failCount = 0;
            Do[
                jExpr = 0; idx = 1;
                Do[
                    alpha = alphas[[ai]];
                    Do[
                        beta = betas[[bi]];
                        coeff = coeffMat[[idx, si+1]];
                        If[coeff != 0,
                            gExp = Table[nuSampleKira[[jj]] - alpha[[jj]], {jj, NE}];
                            term = coeff *
                                Product[nuSampleKira[[jj]]^beta[[jj]], {jj, NE}] *
                                "j" @@ gExp;
                            jExpr = jExpr + term;
                        ];
                        idx++;
                    , {bi, Length[betas]}];
                , {ai, Length[alphas]}];

                reduced = kiraReduce[jExpr /. epsSol];
                If[reduced === 0, passCount++, failCount++];
            , {si, 0, sols-1}];

            kiraTotalPass += passCount; kiraTotalFail += failCount;
            Print["  (lev=", lev, ", deg=", deg,
                ") Kira: ", If[failCount==0, "PASS",
                    "FAIL (" <> ToString[failCount] <> "/" <> ToString[sols] <> ")"],
                "  pass=", passCount, "/", sols];
        , {i, Length[$AllRelations]}];

        Print["\n  KiraVerify: ", kiraTotalPass, "/",
            kiraTotalPass + kiraTotalFail, " PASSED",
            If[kiraTotalFail > 0, " (" <> ToString[kiraTotalFail] <> " failed)", ""]];
    ];
];

If[!TrueQ[$LoadedByBladeFamily], Exit[If[kiraTotalFail > 0, 1, 0]]];
