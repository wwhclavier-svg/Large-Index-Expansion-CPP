(* family_completeness.wl — 特征方程完备性检验(通用版, 由 db313_completeness.wl 推广)
   用法: wolframscript -script tools/family_completeness.wl <family> <order>
   示例: wolframscript -script tools/family_completeness.wl NP222 4

   方法 (Goals.md 完备性标准):
     特征方程 = 替换 v_i -> sector_i*n 后取 n 的首项系数
     (直接从 AllRelations 的 (Alphas,Betas,Coefficients) 构造, 与 characteristicEquation[] 等价)
   对每个含关系的 (lev,deg) 块及累积集合:
     1. 有效性: 特征多项式在该 sector 的 NB=1 regime 点上模 p 求值 (B=Ainv 约定, 应全为0)
     2. 完备性: Singular 求特征理想 dim/vdim (mod p)
        完备 ⟺ dim=0 且 vdim = 该 sector 各 regime 的 Σ NB
   注意: 假设关系针对 top sector (regime sector 分量求和最大者) 重构.
*)

args = Rest[$ScriptCommandLine];
If[Length[args] < 2,
  Print["Usage: wolframscript -script family_completeness.wl <family> <order> [relationsDir]"]; Exit[1]];
{family, order} = {args[[1]], ToExpression[args[[2]]]};

cppRoot = "/home/ykm/work-knowledge/Projects/Large-Index-Expansion-CPP";
relDir = If[Length[args] >= 3, args[[3]], FileNameJoin[{cppRoot, "relations"}]];
levMax = If[Length[args] >= 4, ToExpression[args[[4]]], Infinity];
degMax = If[Length[args] >= 5, ToExpression[args[[5]]], Infinity];
outFile = FileNameJoin[{cppRoot, "results", "intermediate",
  family <> "_completeness_k" <> ToString[order] <> ".json"}];
tmpDir = FileNameJoin[{"/tmp", family <> "_completeness"}];
If[FileType[tmpDir] === None, CreateDirectory[tmpDir]];

Print["Loading relations... ", DateString[]];
Get[FileNameJoin[{relDir, "AllRelations_" <> family <> "_k" <> ToString[order] <> ".m"}]];
Get[FileNameJoin[{relDir, "RelationMeta_" <> family <> ".m"}]];

p = $AllRelations[[1]]["Modulus"];
ne = $AllRelations[[1]]["NE"];

(* top sector = regime sector 分量求和最大 *)
regimes = $RelationMeta["Regimes"];
sectorSums = Total[#["Sector"]] & /@ regimes;
sector = regimes[[First@Ordering[sectorSums, -1]]]["Sector"];
topRegimes = Cases[regimes, r_ /; r["Sector"] == sector];
refVdim = Total[#["NB"] & /@ topRegimes];
validPtsA = Cases[topRegimes, r_ /; r["NB"] == 1 :> r["A"][[;;,1]]];
validPtsAinv = Cases[topRegimes, r_ /; r["NB"] == 1 :> r["Ainv"][[;;,1]]];
Print["family=", family, " k=", order, " ne=", ne, " p=", p];
Print["top sector=", sector, "  regimes: ", Length[topRegimes],
  " (NB=1: ", Length[validPtsA], "),  ref vdim (Σ NB) = ", refVdim];

BV = Array[ToExpression["b" <> ToString[#]] &, ne];

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

clearAssoc[assoc_] := Module[{keys = Keys[assoc], mins},
  If[Length[keys] == 0, Return[assoc]];
  mins = Min /@ Transpose[keys];
  KeyMap[# - mins &, assoc]
];

assocToSingular[assoc_] := StringJoin@Riffle[
  KeyValueMap[
    ToString[#2] <> StringJoin@Table[
      If[#1[[k]] > 0, "*b" <> ToString[k] <> "^" <> ToString[#1[[k]]], ""],
      {k, ne}] &,
    assoc],
  " + "];

evalAssoc[assoc_, pt_] := Mod[Plus @@ KeyValueMap[
  #2 * Times @@ Table[PowerMod[pt[[k]], #1[[k]], p], {k, ne}] &, assoc], p];

singularDim[polys_, tag_] := Module[
  {inF = FileNameJoin[{tmpDir, tag <> ".sing"}],
   outF = FileNameJoin[{tmpDir, tag <> ".out"}], varStr, code, res, dimV, vdimV},
  varStr = StringRiffle["b" <> ToString[#] & /@ Range[ne], ","];
  code = "ring r = " <> ToString[p] <> ", (" <> varStr <> "), dp;\n" <>
    "ideal I = " <> StringRiffle[polys, ",\n"] <> ";\n" <>
    "ideal G = groebner(I);\n" <>
    "printf(\"DIM=%s\", dim(G));\n" <>
    "printf(\"VDIM=%s\", vdim(G));\n" <>
    "ring rl = " <> ToString[p] <> ", (" <> varStr <> ", t), dp;\n" <>
    "ideal I = " <> StringRiffle[polys, ",\n"] <> ";\n" <>
    "ideal J = I + ideal(t*(" <> StringJoin@Riffle["b" <> ToString[#] & /@ Range[ne], "*"] <> ") - 1);\n" <>
    "ideal H = groebner(J);\n" <>
    "printf(\"LVDIM=%s\", vdim(H));\n";
  Export[inF, code, "Text"];
  Run["Singular -q " <> inF <> " > " <> outF <> " 2>&1"];
  res = Import[outF, "Text"];
  dimV = StringCases[res, RegularExpression["(?m)^DIM=(-?\\d+)"] :> "$1"];
  vdimV = StringCases[res, RegularExpression["VDIM=(-?\\d+)"] :> "$1"];
  lvdimV = StringCases[res, RegularExpression["LVDIM=(-?\\d+)"] :> "$1"];
  <|"dim" -> If[dimV =!= {}, ToExpression[dimV[[1]]], "err"],
    "vdim" -> If[vdimV =!= {}, ToExpression[vdimV[[1]]], "err"],
    "lvdim" -> If[lvdimV =!= {}, ToExpression[lvdimV[[1]]], "err"],
    "raw" -> StringTake[res, UpTo[200]]|>
];

blocks = Select[$AllRelations, #["NumRelations"] > 0 && #["Lev"] <= levMax && #["Deg"] <= degMax &];
Print["Blocks with relations (lev<=", levMax, ", deg<=", degMax, "): ",
  {#["Lev"], #["Deg"], #["NumRelations"]} & /@ blocks];

results = <|"family" -> family, "order" -> order, "ne" -> ne, "prime" -> p,
  "topSector" -> sector, "refVdim" -> refVdim|>;
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
  resA = Table[Max[evalAssoc[#, pt] & /@ assocC], {pt, validPtsA}];
  resAinv = Table[Max[evalAssoc[#, pt] & /@ assocC], {pt, validPtsAinv}];
  Print["  validity max residue (B=A):    ", resA];
  Print["  validity max residue (B=Ainv): ", resAinv];
  polys = DeleteCases[assocToSingular /@ assocC, "0" | ""];
  t0 = AbsoluteTime[];
  sd = singularDim[polys, tag];
  Print["  singular dim=", sd["dim"], " vdim=", sd["vdim"], " lvdim=", sd["lvdim"], "  (", Round[AbsoluteTime[]-t0, 0.1], " s)"];
  If[sd["dim"] === "err", Print["  RAW: ", sd["raw"]]];
  results[tag] = <|"lev" -> lev, "deg" -> deg, "numRelations" -> blk["NumRelations"],
    "numCharPolys" -> Length[polys], "validityMaxResidueA" -> resA, "validityMaxResidueAinv" -> resAinv,
    "dim" -> sd["dim"], "vdim" -> sd["vdim"], "lvdim" -> sd["lvdim"],
    "complete" -> (sd["dim"] === 0 && sd["lvdim"] === refVdim && Max[resAinv] === 0)|>;
  cumul = Join[cumul, polys];
  , {blk, blocks}];

Print["\n=== cumulative (all blocks) ===  ", DateString[]];
t0 = AbsoluteTime[];
sd = singularDim[cumul, "cumul"];
Print["  singular dim=", sd["dim"], " vdim=", sd["vdim"], " lvdim=", sd["lvdim"], "  (", Round[AbsoluteTime[]-t0, 0.1], " s)"];
If[sd["dim"] === "err", Print["  RAW: ", sd["raw"]]];
results["cumulative"] = <|"numCharPolys" -> Length[cumul], "dim" -> sd["dim"], "vdim" -> sd["vdim"], "lvdim" -> sd["lvdim"],
  "complete" -> (sd["dim"] === 0 && sd["lvdim"] === refVdim)|>;

Export[outFile, results, "JSON"];
Print["\n-> ", outFile];
Print["Done. ", DateString[]];
