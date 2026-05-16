(* AllRelationsToNuFormat.wl -- Convert AllRelations to Nu-expansion format *)

If[Not@MemberQ[$ContextPath, "AllRelationsToNuFormat`"],
    PrependTo[$ContextPath, "AllRelationsToNuFormat`"]
];

BeginPackage["AllRelationsToNuFormat`"]

AllRelationsToNuFormat::usage = "AllRelationsToNuFormat[$AllRelations] or AllRelationsToNuFormat[$AllRelations, Gsymbol, vlist] converts C++ AllRelations to NuFormat Association.\n\
Returns <|{lev,deg} -> Association|>.\n\
Optional: Gsymbol (default G), vlist (default {nu1,...,nuNE})."

AllRelationsToNuFormatFile::usage = "AllRelationsToNuFormatFile[filepath] loads AllRelations_*.m and converts it."

Begin["`Private`"]

buildNuRelation[entry_, solCol_, nuSyms_, Gsymbol_] := Module[
    {alphas, betas, coeffMat, nAlpha, nBeta, idx, alpha, beta,
     coeff, betaPoly, terms},
    alphas = entry["Alphas"]; betas = entry["Betas"];
    coeffMat = entry["Coefficients"];
    nAlpha = Length[alphas]; nBeta = Length[betas]; idx = 1; terms = {};
    Do[alpha = alphas[[ai]];
       Do[beta = betas[[bi]];
          coeff = coeffMat[[idx, solCol + 1]];
          If[coeff != 0,
             betaPoly = If[AllTrue[beta, # == 0 &], 1,
                 Times @@ Thread[Power[nuSyms, beta]]];
             AppendTo[terms, coeff * betaPoly * Gsymbol @@ Thread[nuSyms - alpha]];
          ];
          idx++;, {bi, nBeta}];, {ai, nAlpha}];
    Which[Length[terms] == 0, 0, Length[terms] == 1, First[terms], True, Plus @@ terms]
];

AllRelationsToNuFormat[allRelations_List, Gsymbol_Symbol : G, vlist_ : Automatic] := Module[
    {NE, nuSyms, modulus, result, entry, lev, deg, nSols, key,
     relParticular, basisRels},
    NE = allRelations[[1]]["NE"];
    nuSyms = If[vlist === Automatic, Table[Symbol["nu" <> ToString[i]], {i, NE}], vlist];
    modulus = allRelations[[1]]["Modulus"];
    result = Association[];
    Do[
        entry = allRelations[[i]];
        {lev, deg, nSols} = entry /@ {"Lev", "Deg", "NumSolutions"};
        key = {lev, deg};
        relParticular = buildNuRelation[entry, 0, nuSyms, Gsymbol];
        basisRels = Table[buildNuRelation[entry, col, nuSyms, Gsymbol], {col, 1, nSols}];
        result[key] = <|
            "Relations" -> basisRels,
            "Particular" -> relParticular, "Basis" -> basisRels,
            "NE" -> NE, "Modulus" -> modulus,
            "Alphas" -> entry["Alphas"], "Betas" -> entry["Betas"]
        |>;
    , {i, Length[allRelations]}];
    result
];

AllRelationsToNuFormatFile[filepath_String] := Module[{},
    Get[filepath];
    If[!ValueQ[Global`$AllRelations],
        Print["ERROR: $AllRelations not defined"]; Return[$Failed]];
    AllRelationsToNuFormat[Global`$AllRelations]
]

End[]

EndPackage[]

(* -- CLI -- *)
If[Length[$ScriptCommandLine] >= 2 && $InputFileName === $ScriptCommandLine[[1]],
    Module[{famname, projectRoot, relationsDir, files, pickK, f, result, out, content},
        famname = $ScriptCommandLine[[2]];
        projectRoot = DirectoryName[ExpandFileName[$InputFileName]] <> "/../../..";
        relationsDir = projectRoot <> "/relations";
        files = Select[FileNames["AllRelations_" <> famname <> "_k*.m", relationsDir],
            StringMatchQ[FileNameTake[#], "AllRelations_" <> famname <> "_k" ~~ DigitCharacter .. ~~ ".m"] &];
        If[files === {}, Print["ERROR: No AllRelations_", famname, "_k*.m"]; Exit[1]];
        pickK[f_] := ToExpression[StringCases[FileNameTake[f], "k" ~~ d : DigitCharacter .. ~~ ".m" :> d][[1]]];
        f = Last[SortBy[files, pickK]];
        Print["Loading: ", FileNameTake[f]];
        result = AllRelationsToNuFormatFile[f];
        If[result === $Failed, Exit[1]];
        out = FileNameJoin[{projectRoot, "NuFormat_" <> famname <> ".m"}];
        Print["Writing: ", FileNameTake[out]];
        Save[out, result];
        content = Import[out, "Text"];
        Export[out, StringReplace[content, "result = " -> "$NuFormat = "], "Text"];
        Print["Done: ", out];
        Exit[0];
    ]
];
