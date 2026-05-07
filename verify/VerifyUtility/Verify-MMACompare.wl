(* ::Package:: *)
(* Verify-MMACompare.wl — MMA relation coefficient comparison *)
(* Usage standalone: wolframscript -file Verify-MMACompare.wl <famname> *)
(* Load in-process:  $LoadedByBladeFamily=True; Get["Verify-MMACompare.wl"] *)

If[!TrueQ[$LoadedByBladeFamily],
    current = If[$FrontEnd === Null, $InputFileName, NotebookFileName[]] // DirectoryName;
    $ProjectRoot = FileNameJoin[{current, ".."}];
    Get[FileNameJoin[{$ProjectRoot, "verify/FamilyDatabase/FamilyDatabase.wl"}]];

    If[Length[$ScriptCommandLine] >= 2, famname = $ScriptCommandLine[[2]], famname = "bub00"];
    If[!KeyExistsQ[$FamilyDatabase, famname],
        Print["ERROR: Unknown family '", famname, "'"]; Exit[1]];

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
];

Print["\n=== MODULE 4: MMACompare ==============================="];

$MMAFile = FileNameJoin[{$ProjectRoot, "MMARelations_" <> famname <> ".m"}];
If[!FileExistsQ[$MMAFile],
    Print["  SKIP: ", FileNameTake[$MMAFile], " not found"];
,
    Get[$MMAFile];
    If[!ValueQ[$MMARelations],
        Print["  SKIP: $MMARelations not defined"];
        Print["  Expected: $MMARelations = <| {lev,deg} -> <|\"exprs\"->{...}|> |>"];
    ,
        mmaTotalPass = 0; mmaTotalFail = 0;
        Do[
            entry = $AllRelations[[i]];
            lev = entry["Lev"]; deg = entry["Deg"]; sols = entry["NumSolutions"];
            If[sols == 0 || !KeyExistsQ[$MMARelations, {lev, deg}], Continue[]];

            mmaExprs = $MMARelations[{lev, deg}]["exprs"];
            Do[
                Print["  (lev=", lev, ", deg=", deg, ") Sol ", si,
                    " MMA: not implemented"];
            , {si, 1, Min[Length[mmaExprs], sols]}];
        , {i, Length[$AllRelations]}];

        Print["\n  MMACompare: see above (custom implementation needed)"];
    ];
];

If[!TrueQ[$LoadedByBladeFamily], Exit[0]];
