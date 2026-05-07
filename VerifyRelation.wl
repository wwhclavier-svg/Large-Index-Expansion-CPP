(* ::Package:: *)
(* VerifyRelation.wl — Unified LIE Relation Verification *)
(* Usage: wolframscript -file VerifyRelation.wl <famname> [nuSample] [skip M1,M2] *)
(*
    Four modules, run in-process (shared data loaded once):
      Module 1: BladeVerify   — Blade IBP reduction at test nu points
      Module 2: SeriesVerify   — Asymptotic expansion order-by-order substitution
      Module 3: KiraVerify     — Kira IBP rule cross-validation
      Module 4: MMACompare     — MMA relation coefficient comparison
*)

(* ====================== Path Initialization ====================== *)
current = If[$FrontEnd === Null, $InputFileName, NotebookFileName[]] // DirectoryName;
$ProjectRoot = current;
$BladePath = "/home/ykm/blade-workspace";
$VerifyDir = FileNameJoin[{$ProjectRoot, "verify", "VerifyUtility"}];

(* Load Blade first to avoid symbol shadowing with FamilyDatabase *)
Get[FileNameJoin[{$BladePath, "Blade.wl"}]];

Get[FileNameJoin[{$ProjectRoot, "verify/FamilyDatabase/FamilyDatabase.wl"}]];

(* ====================== CLI Argument Parsing ====================== *)
If[Length[$ScriptCommandLine] >= 2,
    famname = $ScriptCommandLine[[2]],
    famname = "bub00"];
nuSampleOverride = If[Length[$ScriptCommandLine] >= 3 &&
    StringMatchQ[$ScriptCommandLine[[3]], "{" ~~ __],
    ToExpression[$ScriptCommandLine[[3]]],
    Automatic];
skipModules = {};
Do[
    If[$ScriptCommandLine[[i]] == "skip" && i < Length[$ScriptCommandLine],
        skipModules = StringSplit[$ScriptCommandLine[[i+1]], ","]],
    {i, Length[$ScriptCommandLine]}];

If[!KeyExistsQ[$FamilyDatabase, famname],
    Print["ERROR: Unknown family '", famname, "'"];
    Print["Available: ", Keys[$FamilyDatabase]];
    Exit[1]];

(* ====================== Family Configuration ====================== *)
config = $FamilyDatabase[famname];
bladeName = StringReplace[famname, "-" -> ""];
family = Symbol[bladeName];

dimension = 4 - 2*eps;
loop = config["LoopMomenta"];
leg = config["ExternalMomenta"];
conservation = {};
replacement = config["KinematicRules"];
propagator = config["Propagators"];
topsector = config["TopSector"];
numeric = DeleteCases[config["Numeric"], "d" -> _];
dValue = "d" /. config["Numeric"];
modulus = config["Modulus"];

Print["====================================================="];
Print["  VerifyRelation: ", famname];
Print["  ", config["Description"]];
Print["  TopSector: ", topsector, "  Props: ", Length[propagator],
      "  Loops: ", Length[loop], "  d: ", dValue];
If[skipModules =!= {}, Print["  Skipping: ", StringRiffle[skipModules, ", "]]];
Print["====================================================="];

(* Blade family definition — only if not skipped *)
If[!MemberQ[skipModules, "BladeVerify"],
    BLSetReducerOptions["MaxIncrement" -> 5, "BlackBoxRank" -> 3, "BlackBoxDot" -> 1];
    BLFamilyDefine[family, dimension, propagator, loop, leg,
        conservation, replacement, topsector, numeric];
    epsSol = Solve[dimension - dValue == 0, eps][[1]];
];

(* ====================== Load C++ AllRelations ====================== *)
allRelFiles = FileNames["AllRelations_" <> famname <> "_k*.m", $ProjectRoot];
If[allRelFiles === {},
    Print["ERROR: No AllRelations_", famname, "_k*.m found"]; Exit[1]];
allRelFilesNum = Select[allRelFiles,
    StringMatchQ[FileNameTake[#],
        "AllRelations_" <> famname <> "_k" ~~ DigitCharacter .. ~~ ".m"] &];
If[allRelFilesNum === {},
    Print["ERROR: No AllRelations_", famname, "_k<N>.m found"]; Exit[1]];
extractK[fname_String] := ToExpression[
    StringCases[fname, "AllRelations_" <> famname <> "_k" ~~ k : DigitCharacter .. ~~ ".m" :> k][[1]]];
allRelFile = Last[SortBy[allRelFilesNum, extractK]];
Print["\n  AllRelations: ", FileNameTake[allRelFile]];
Get[allRelFile];
If[!ValueQ[$AllRelations],
    Print["ERROR: $AllRelations not defined"]; Exit[1]];

NE = $AllRelations[[1]]["NE"];
mkNuSample[] := If[nuSampleOverride =!= Automatic, nuSampleOverride,
    If[NE <= 2, {2, 5}, If[NE <= 5, {1, 2, 1, 2, 1}, Table[RandomInteger[{1, 15}], {NE}]]]];

(* Config summary *)
nConfigs = Length[$AllRelations];
nNonTrivial = Length[Select[$AllRelations, #["NumSolutions"] > 0 &]];
Print["  NE: ", NE, "  Configs: ", nConfigs, " (", nNonTrivial, " non-trivial)"];
Do[
    Print["  (lev=", StringPadRight[ToString[#["Lev"]], 3],
        " deg=", StringPadRight[ToString[#["Deg"]], 3],
        ")  sols=", StringPadRight[ToString[#["NumSolutions"]], 4],
        " |alpha|=", StringPadRight[ToString[Length[#["Alphas"]]], 5],
        " |beta|=", StringPadRight[ToString[Length[#["Betas"]]], 5],
        " stable=", #["StableOrder"]] & /@ $AllRelations;
];

(* ====================== Load Modules In-Process ====================== *)
$LoadedByBladeFamily = True;

If[!MemberQ[skipModules, "BladeVerify"],
    Get[FileNameJoin[{$VerifyDir, "Verify-Blade.wl"}]]];

If[!MemberQ[skipModules, "SeriesVerify"],
    Get[FileNameJoin[{$VerifyDir, "Verify-Series.wl"}]]];

If[!MemberQ[skipModules, "KiraVerify"],
    Get[FileNameJoin[{$VerifyDir, "Verify-Kira.wl"}]]];

If[!MemberQ[skipModules, "MMACompare"],
    Get[FileNameJoin[{$VerifyDir, "Verify-MMACompare.wl"}]]];

(* ====================== Summary ====================== *)
Print["\n\n====================================================="];
Print["  SUMMARY: ", famname];
Print["====================================================="];
If[ValueQ[bladeTotalPass],
    Print["  BladeVerify:   ", bladeTotalPass, "/", bladeTotalPass + bladeTotalFail,
        If[bladeTotalFail > 0, " (" <> ToString[bladeTotalFail] <> " failed)", " PASS"]]];
If[ValueQ[seriesTotalPass],
    Print["  SeriesVerify:  ",
        If[seriesTotalPass + seriesTotalFail > 0,
            ToString[seriesTotalPass] <> "/" <> ToString[seriesTotalPass + seriesTotalFail] <>
            If[seriesTotalFail > 0,
                " (" <> ToString[seriesTotalFail] <> " failed)", " PASS"], "N/A"]]];
If[ValueQ[kiraTotalPass],
    Print["  KiraVerify:    ",
        If[kiraTotalPass + kiraTotalFail > 0,
            ToString[kiraTotalPass] <> "/" <> ToString[kiraTotalPass + kiraTotalFail] <>
            If[kiraTotalFail > 0,
                " (" <> ToString[kiraTotalFail] <> " failed)", " PASS"], "N/A"]]];
Print["====================================================="];
