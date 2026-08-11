(* DB313 k3: characteristic-equation completeness check
   直接构造法: 特征方程 = 替换 v_i -> sector_i*n 后取 n 的首项系数
   (与 characteristicEquation[] 数学等价, 但不构造巨型符号表达式)
   对每个含关系的 (lev,deg) 块及累积集合:
     1. 有效性: 特征多项式在 RelationMeta 的 7 个 NB=1 regime 点上模 p 求值(应全为0)
     2. 完备性: Singular 求特征理想的 dim / vdim (模 p)
   参考: docs/LIE_Benchmark_Record.md 第5节, Goals.md 完备性标准
*)

p = 179424673;
ne = 9;
sector = {1,1,1,1,1,1,1,0,0};
cppRoot = "/home/ykm/work-knowledge/Projects/Large-Index-Expansion-CPP";
outFile = FileNameJoin[{cppRoot, "results", "intermediate", "db313_completeness_k3.json"}];
tmpDir = "/tmp/db313_completeness";
If[FileType[tmpDir] === None, CreateDirectory[tmpDir]];

Print["Loading relations... ", DateString[]];
Get[FileNameJoin[{cppRoot, "relations", "AllRelations_DB313_k3.m"}]];
Get[FileNameJoin[{cppRoot, "relations", "RelationMeta_DB313.m"}]];
Print["Loaded. Memory: ", MemoryInUse[]/2.^20, " MB"];

BV = Array[ToExpression["b" <> ToString[#]] &, ne];

(* ---- 直接构造特征多项式: 返回 {<expvec -> coeff, ...>} (mod p) ---- *)
buildCharAssoc[entry_] := Module[
  {alphas = entry["Alphas"], betas = entry["Betas"], coefs = entry["Coefficients"],
   nA, nB, nR, w, wmax, sel},
  nA = Length[alphas]; nB = Length[betas]; nR = entry["NumRelations"];
  If[Dimensions[coefs] =!= {nA*nB, nR+1},
    Print["WARN dims mismatch: ", Dimensions[coefs], " vs ", {nA*nB, nR+1}]];
  w = (sector . #) & /@ betas;
  wmax = Max[w];
  sel = Flatten@Position[w, wmax];
  Table[
    Merge[
      DeleteCases[
        Flatten@Table[
          If[coefs[[(i-1)*nB+j, r+1]] =!= 0, alphas[[i]] -> coefs[[(i-1)*nB+j, r+1]], Nothing],
          {j, sel}, {i, nA}],
        _ -> 0],
      Total],
    {r, nR}]
];

(* 清分母: 每个多项式减去逐变量最小指数 *)
clearAssoc[assoc_] := Module[{keys = Keys[assoc], mins},
  If[Length[keys] == 0, Return[assoc]];
  mins = Min /@ Transpose[keys];
  KeyMap[# - mins &, assoc]
];

assocToPoly[assoc_] := Plus @@ KeyValueMap[
  #2 * Times @@ (BV^#1) &, assoc];

assocToSingular[assoc_] := StringJoin@Riffle[
  KeyValueMap[
    ToString[#2] <> StringJoin@Table[
      If[#1[[k]] > 0, "*b" <> ToString[k] <> "^" <> ToString[#1[[k]]], ""],
      {k, ne}] &,
    assoc],
  " + "];

(* ---- 有效性检查点: NB=1 regimes; B=A 或 B=Ainv 两个约定都试 ---- *)
validPtsA = Cases[$RelationMeta["Regimes"], r_ /; r["NB"] == 1 :> r["A"][[;;,1]]];
validPtsAinv = Cases[$RelationMeta["Regimes"], r_ /; r["NB"] == 1 :> r["Ainv"][[;;,1]]];
Print["Validity points (NB=1 regimes): ", Length[validPtsA]];

evalAssoc[assoc_, pt_] := Mod[Plus @@ KeyValueMap[
  #2 * Times @@ Table[PowerMod[pt[[k]], #1[[k]], p], {k, ne}] &, assoc], p];

(* ---- Singular dim/vdim ---- *)
singularDim[polys_, tag_] := Module[
  {inF = FileNameJoin[{tmpDir, tag <> ".sing"}],
   outF = FileNameJoin[{tmpDir, tag <> ".out"}], varStr, code, res},
  varStr = StringRiffle["b" <> ToString[#] & /@ Range[ne], ","];
  code = "ring r = " <> ToString[p] <> ", (" <> varStr <> "), dp;\n" <>
    "ideal I = " <> StringRiffle[polys, ",\n"] <> ";\n" <>
    "ideal G = groebner(I);\n" <>
    "printf(\"DIM=%s\", dim(G));\n" <>
    "printf(\"VDIM=%s\", vdim(G));\n";
  Export[inF, code, "Text"];
  Run["Singular -q " <> inF <> " > " <> outF <> " 2>&1"];
  res = Import[outF, "Text"];
  dimV = StringCases[res, RegularExpression["(?m)^DIM=(-?\\d+)"] :> "$1"];
  vdimV = StringCases[res, RegularExpression["VDIM=(-?\\d+)"] :> "$1"];
  <|"dim" -> If[dimV =!= {}, ToExpression[dimV[[1]]], "err"],
    "vdim" -> If[vdimV =!= {}, ToExpression[vdimV[[1]]], "err"],
    "raw" -> StringTake[res, UpTo[200]]|>
];

(* ---- 主流程 ---- *)
blocks = Select[$AllRelations, #["NumRelations"] > 0 &];
Print["Blocks with relations: ", {#["Lev"], #["Deg"], #["NumRelations"]} & /@ blocks];

results = <||>;
cumul = {};
Do[
  lev = blk["Lev"]; deg = blk["Deg"];
  tag = "l" <> ToString[lev] <> "d" <> ToString[deg];
  Print["\n=== (lev=", lev, ", deg=", deg, ") ===  ", DateString[]];
  t0 = AbsoluteTime[];
  assoc = buildCharAssoc[blk];
  assoc = DeleteCases[assoc, <||>];
  Print["  char polys: ", Length[assoc], "  (build ", Round[AbsoluteTime[]-t0, 0.1], " s)"];
  assocC = clearAssoc /@ assoc;
  (* 有效性: B=A 与 B=Ainv 两种约定 *)
  resA = Table[Max[evalAssoc[#, pt] & /@ assocC], {pt, validPtsA}];
  resAinv = Table[Max[evalAssoc[#, pt] & /@ assocC], {pt, validPtsAinv}];
  Print["  validity max residue (B=A):    ", resA];
  Print["  validity max residue (B=Ainv): ", resAinv];
  (* 完备性 *)
  polys = DeleteCases[assocToSingular /@ assocC, "0" | ""];
  t0 = AbsoluteTime[];
  sd = singularDim[polys, tag];
  Print["  singular dim=", sd["dim"], " vdim=", sd["vdim"], "  (", Round[AbsoluteTime[]-t0, 0.1], " s)"];
  If[sd["dim"] === "err", Print["  RAW: ", sd["raw"]]];
  results[tag] = <|"lev" -> lev, "deg" -> deg, "numRelations" -> blk["NumRelations"],
    "numCharPolys" -> Length[polys], "validityMaxResidueA" -> resA, "validityMaxResidueAinv" -> resAinv,
    "dim" -> sd["dim"], "vdim" -> sd["vdim"]|>;
  cumul = Join[cumul, polys];
  , {blk, blocks}];

(* 累积集合 *)
Print["\n=== cumulative (all blocks) ===  ", DateString[]];
t0 = AbsoluteTime[];
sd = singularDim[cumul, "cumul"];
Print["  singular dim=", sd["dim"], " vdim=", sd["vdim"], "  (", Round[AbsoluteTime[]-t0, 0.1], " s)"];
If[sd["dim"] === "err", Print["  RAW: ", sd["raw"]]];
results["cumulative"] = <|"numCharPolys" -> Length[cumul], "dim" -> sd["dim"], "vdim" -> sd["vdim"]|>;

Export[outFile, results, "JSON"];
Print["\n-> ", outFile];
Print["Done. ", DateString[]];
