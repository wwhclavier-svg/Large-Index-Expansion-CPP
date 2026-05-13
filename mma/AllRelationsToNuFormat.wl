(* ::Package:: *)
(* AllRelationsToNuFormat.wl — Convert C++ AllRelations to ν-expansion format *)
(*
   Input:  AllRelations_<fam>_k*.m  (C++ exported Association list)
   Output: NuFormat_<fam>.m          ($NuFormat Association, each relation as
                                     Σ b(α,β) · ν^β · G[ν-α]  form)

   Usage:
     wolframscript -file AllRelationsToNuFormat.wl <famname>

   Output format:
     $NuFormat = <|
       {lev, deg} -> <|
         "Relations" -> { relParticular, relBasis1, relBasis2, ... },
         "Particular" -> relParticular,    (* particular solution *)
         "Basis" -> { relBasis1, ... },    (* null space basis vectors *)
         "NE" -> 2,
         "Modulus" -> 179424673,
         "Alphas" -> {{0,0}, ...},
         "Betas"  -> {{0,0}, ...}
       |>,
       ...
     |>

   Each relation is a sum of terms:
     ∑_{α∈Alphas, β∈Betas} b(α,β) · (ν_1^{β_1} ... ν_NE^{β_NE}) · G[ν-α]

   where G[ν-α] = G[ν1-α1, ν2-α2, ...].

   For NumSolutions == 0: only Particular exists, Basis = {}.
   For NumSolutions  > 0: Particular uses arbitrary values for free vars (=0),
                          Basis[k] is the k-th null space direction.
*)

(* ---- Parse CLI ---- *)
If[Length[$ScriptCommandLine] < 2,
    Print["Usage: wolframscript -file AllRelationsToNuFormat.wl <famname>"];
    Print["Example: wolframscript -file AllRelationsToNuFormat.wl bub00"];
    Exit[1];
];

famname = $ScriptCommandLine[[2]];
projectRoot = DirectoryName[$InputFileName] <> "/..";
relationsDir = projectRoot <> "/relations";

(* ---- Load AllRelations ---- *)
allRelFiles = FileNames["AllRelations_" <> famname <> "_k*.m", relationsDir];
If[allRelFiles === {},
    Print["ERROR: No AllRelations_", famname, "_k*.m found in ", relationsDir];
    Exit[1];
];

(* Pick largest k *)
allRelFilesNum = Select[allRelFiles,
    StringMatchQ[FileNameTake[#],
        "AllRelations_" <> famname <> "_k" ~~ DigitCharacter .. ~~ ".m"] &];
extractK[fname_String] := ToExpression[
    StringCases[fname, "AllRelations_" <> famname <> "_k" ~~ k : DigitCharacter .. ~~ ".m" :> k][[1]]];
allRelFile = Last[SortBy[allRelFilesNum, extractK]];

Print["Loading: ", FileNameTake[allRelFile]];
Get[allRelFile];

If[!ValueQ[$AllRelations],
    Print["ERROR: $AllRelations not defined after loading"];
    Exit[1];
];

Print["  Configs loaded: ", Length[$AllRelations]];
Print["  NE = ", $AllRelations[[1]]["NE"]];
Print["  Modulus = ", $AllRelations[[1]]["Modulus"]];

(* ---- Build ν symbols (plain ASCII for MMA CLI compatibility) ---- *)
NE = $AllRelations[[1]]["NE"];
nuSyms = Table[Symbol["nu" <> ToString[i]], {i, NE}];
modulus = $AllRelations[[1]]["Modulus"];

(* ---- Build each relation ---- *)
buildNuRelation[entry_, solCol_] := Module[
    {alphas, betas, coeffMat, nAlpha, nBeta, idx, terms, alpha, beta,
     coeff, betaPoly, gArgs, term},

    alphas = entry["Alphas"];
    betas  = entry["Betas"];
    coeffMat = entry["Coefficients"];
    nAlpha = Length[alphas];
    nBeta  = Length[betas];
    idx = 1;
    terms = {};

    Do[
        alpha = alphas[[ai]];
        Do[
            beta = betas[[bi]];
            coeff = coeffMat[[idx, solCol + 1]];  (* 1-indexed: col 1 = particular *)
            If[coeff != 0,
                (* ν^β = ν_1^{β_1} * ν_2^{β_2} * ... *)
                betaPoly = Times @@ Thread[Power[nuSyms, beta]];
                If[betaPoly === 1, betaPoly = 1];

                (* G[ν-α] *)
                gArgs = Thread[nuSyms - alpha];
                term = coeff * betaPoly * G @@ gArgs;
                AppendTo[terms, term];
            ];
            idx++;
        , {bi, nBeta}];
    , {ai, nAlpha}];

    If[Length[terms] == 0, Return[0]];
    If[Length[terms] == 1, Return[First[terms]]];
    Plus @@ terms
];

(* ---- Process all entries ---- *)
result = Association[];

Do[
    entry = $AllRelations[[i]];
    lev = entry["Lev"];
    deg = entry["Deg"];
    nSols = entry["NumSolutions"];

    key = {lev, deg};

    (* Particular solution (col 0): unique relation if nSols==0, arbitrary if nSols>0 *)
    relParticular = buildNuRelation[entry, 0];

    (* Null space basis vectors *)
    basisRels = Table[buildNuRelation[entry, col], {col, 1, nSols}];

    allRels = Join[{relParticular}, basisRels];

    result[key] = <|
        "Relations" -> allRels,
        "Particular" -> relParticular,
        "Basis" -> basisRels,
        "NE" -> NE,
        "Modulus" -> modulus,
        "Alphas" -> entry["Alphas"],
        "Betas"  -> entry["Betas"]
    |>;

    nTermsParticular = Which[
        relParticular === 0, 0,
        Head[relParticular] === Plus, Length[relParticular],
        True, 1];
    nTermsBasis = Total[Which[
        # === 0, 0,
        Head[#] === Plus, Length[#],
        True, 1] & /@ basisRels];
    Print["  (lev=", lev, ", deg=", deg, ")  nSols=", nSols,
    "  terms(part)=", nTermsParticular,
    "  terms(basis)=", nTermsBasis];
, {i, Length[$AllRelations]}];

(* ---- Write output ---- *)
outputFile = FileNameJoin[{projectRoot, "NuFormat_" <> famname <> ".m"}];
Print["\nWriting: ", FileNameTake[outputFile]];
Print["  Total configs: ", Length[Keys[result]]];
Print["  Output variable: $NuFormat"];

(* Write as .m file *)
Save[outputFile, result];
(* Also export the variable name so it loads as $NuFormat *)
content = Import[outputFile, "Text"];
(* Replace "result = " with "$NuFormat = " *)
content = StringReplace[content, "result = " -> "$NuFormat = "];
Export[outputFile, content, "Text"];

Print["\nDone. Output: ", outputFile];
