(* ::Package:: *)

(* blade_test_fixed.m \[LongDash] \:4fee\:590d\:7248 *)
(* \:8fd0\:884c: cd build && wolframscript -file blade_test_fixed.m *)

$ProjectRoot = DirectoryName[$InputFileName] <> "/..";
current = If[$FrontEnd===Null,$InputFileName,NotebookFileName[]]//DirectoryName;

$BladePath = "/home/ykm/blade-workspace";
Get[FileNameJoin[{$BladePath, "Blade.wl"}]];

Get[FileNameJoin[{current, "../verify/FamilyDatabase/FamilyDatabase.wl"}]];


famname = "bub00";
family = Symbol[StringReplace[famname, "-" -> ""]];  (* bub00 *)

(* \:4ece FamilyDatabase \:52a0\:8f7d\:914d\:7f6e *)
dimension = 4 - 2*eps;
loop = $FamilyDatabase[famname]["LoopMomenta"];
leg = $FamilyDatabase[famname]["ExternalMomenta"];
conservation = {};
replacement = $FamilyDatabase[famname]["KinematicRules"];
propagator = $FamilyDatabase[famname]["Propagators"];
topsector = $FamilyDatabase[famname]["TopSector"];
numeric = $FamilyDatabase[famname]["Numeric"];
modulus = $FamilyDatabase[famname]["Modulus"];

Print["numeric: ", numeric, "  modulus: ", modulus];
NE=Length@propagator;Alist=Array[A,Length@propagator];vlist=Array["v"<>ToString[#]&,Length@propagator];

(* \:5b9a\:4e49 Blade \:65cf *)
BLFamilyDefine[family, dimension, propagator, loop, leg,
  conservation, replacement, topsector, numeric];



(* \:52a0\:8f7d C++ \:5173\:7cfb *)
$RelationFile="AllRelations_"<>famname<>"_k6.m";
allRelFile = FileNameJoin[{current, "..",$RelationFile}];
Get[allRelFile];

Print[levdeg={#["Lev"], #["Deg"]} & /@ $AllRelations];
Print[Length /@ {#["Alphas"], #["Betas"]} & /@ $AllRelations];
Print[matdim=Dimensions[#["Coefficients"]] & /@ $AllRelations];
Print[MatrixForm@SparseArray[Thread[Map[#+{0,1}&,levdeg]->matdim[[;;,2]]-1]//DeleteCases[#,a_/;a[[1,1]]==0]&]];

(* \:6784\:5efa relcpp \:5173\:7cfb\:8868 *)
relcpp = Association@Table[
    Block[{lev, deg, alphas, betas, rel},
      {lev, deg} = $AllRelations[[i]] // {#["Lev"], #["Deg"]} &;
      {alphas, betas} = $AllRelations[[i]] // {#["Alphas"], #["Betas"]} &;
      rel = Sum[
        Sum[$AllRelations[[i]]["Coefficients"][[(j-1)*Length[betas] + k]] *
            Times @@ Thread[Power[vlist, betas[[k]]]],
          {k, Length[betas]}] "j" @@ (vlist - alphas[[j]]),
        {j, Length[alphas]}];
      {lev, deg} -> Rest[rel]
    ],
    {i, Length[$AllRelations]}];



nuSample = {2, 5};

(* ---- BL \:7ea6\:5316 ---- *)
Print["j-target: "];
jvar = Cases[relcpp, a_ /; Head[a] === "j", Infinity] // DeleteDuplicates;
jTarget = jvar /. Thread[vlist -> nuSample];
Print[jTarget];

Print["BL-Reduction: "];
res = BLReduce[jTarget /. "j"[a__] :> BL[family, {a}]];

bladeTable = Thread[jTarget -> (res /. BL[family, a__] :> "j" @@ a)];

Print["BL-Reduction Rule: "];
bladeRule = (bladeTable /.
    Solve[dimension - ("d" /. $FamilyDatabase[famname]["Numeric"]) == 0, eps][[1]]) //
  Map[(#[[1]] -> PolynomialMod[#[[2]], modulus]) &, #] &;

(* ---- C++ \:5173\:7cfb\:9a8c\:8bc1 ---- *)
Print["C++ Verification: "];
relcppVeri = (Values[relcpp] /. Thread[vlist -> nuSample] /. bladeRule) //
    PolynomialMod[#, modulus] &;
Print[relcppVeri // MatrixForm];



(* ---- relmma \:5b9a\:4e49\:ff08\:9700\:8981\:624b\:52a8\:586b\:5199 MMA \:7ed3\:679c\:ff09---- *)
(* TODO: \:586b\:5165 MMA \:8ba1\:7b97\:7684\:5173\:7cfb\:8868\:ff0c\:683c\:5f0f\:4e3a 2\[Times]N \:77e9\:9635 *)
(* \:793a\:4f8b\:ff1abub00 (lev=1, deg=1) \:7684\:5173\:7cfb *)
relmma = {
  {(59808225+v1) "j"[-1+v1, v2] + (59808226+179424670 v1+179424671 v2) "j"[v1, -1+v2] + (3+3 v1+179424667 v2) "j"[v1, v2],
   (179424672+v2) "j"[-1+v1, v2] + (59808223+2 v1+v2) "j"[v1, -1+v2] + (179424670+3 v2) "j"[v1, v2]}
};



(* ---- MMA \:5173\:7cfb\:9a8c\:8bc1\:ff08\:5982\:679c relmma \:5df2\:5b9a\:4e49\:ff09---- *)
If[ValueQ[relmma],
  (* \:5c06 symbol v1,v2 \:8f6c\:4e3a string "v1","v2"\:ff0c\:4e0e relcpp \:7ea6\:5b9a\:4e00\:81f4 *)
  relmmaStr = relmma /. ((ToExpression[#]->#)&/@vlist);
  Print["MMA Verification: "];
  relmmaVeri = (relmmaStr /. Thread[vlist -> nuSample] /. bladeRule) //
      PolynomialMod[#, modulus] &;
  Print[relmmaVeri // MatrixForm];
];

(* ---- Series \:9a8c\:8bc1\:ff08\:9700\:8981 ExpansionMMA_bub00.m\:ff09---- *)
$ExpansionFile="ExpansionMMA_"<>famname<>".m";
If[FileExistsQ[FileNameJoin[{$ProjectRoot, $ExpansionFile}]],
  Get[FileNameJoin[{$ProjectRoot, $ExpansionFile}]];
  order = 4;
  hlist = $ExpansionResults[[1, 1]]["Solutions"][[1]]["H"] //
      Sum[1/n^k * (#[[k+1]] /. ((ToExpression[#]->#)&/@vlist)), {k, 0, Length[#]-1}] &;
  secpos = {1, 1};
  vlist = {"v1", "v2"};
  Alist = Array[A, 2];
  varRule = {A[2] -> 59808223, A[1] -> 59808223};

  seriesVerify[relations_, secpos_, degree_, order_] := Module[{rel},
    rel = relations /. {"j"[a__] :> jshift @@ (vlist - {a})};
    rel = (rel / n^degree) /. Thread[vlist -> vlist + secpos*n];
    rel = rel /. {jshift[a__] :> (PolynomialMod[#, modulus] &@
        Times @@ Thread[Power[Alist /. varRule, -{a}]] *
        (hlist /. Thread[vlist -> vlist - {a}]))};
    rel = rel // CoefficientList[#, 1/n, order+1] & // PolynomialMod[#, modulus] &];

  Print["Series Verification: "];
  seriesVerify[relmma, secpos, 1, order] // Print;
,
  Print["SKIP: "<>$ExpansionFile<>" not found \[LongDash] series check skipped."];
];

Print["Done."];
