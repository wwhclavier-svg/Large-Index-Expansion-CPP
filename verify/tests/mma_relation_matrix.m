(* mma_relation_matrix.m — 生成 (lev,deg) 关系数量矩阵 *)
SetDirectory["/root/Large-Index-Expansion-MMA-Mini"];

Get["LIEReconstruct.wl"];
Get["LIEFamilyDefine.wl"];
Get["LIERegions.wl"];
Get["LIEExpand.wl"];
Get["LIEWorkflow.wl"];

modulus = Prime[10000000]; (* 179424673 *)

families = <|
  "bub00" -> <|
    "props" -> {k1^2, -(k1 - p1)^2},
    "numeric" -> {s -> 3, msq -> 0, "d" -> 1/3}
  |>,
  "bub10" -> <|
    "props" -> {-k1^2 + msq, -(k1 - p1)^2},
    "numeric" -> {s -> 3, msq -> 1, "d" -> 1/3}
  |>,
  "bub11" -> <|
    "props" -> {-k1^2 + msq, -(k1 - p1)^2 + msq},
    "numeric" -> {s -> 3, msq -> 1, "d" -> 1/3}
  |>
|>;

Do[
  fam = famName;
  cfg = families[fam];
  Print["========================================"];
  Print["Family: ", fam];
  Print["========================================"];

  props = cfg["props"] /. cfg["numeric"];
  kinRules = {p1^2 -> s} /. cfg["numeric"];

  data = LIEDefineFamily[props, {k1}, {p1}, kinRules, {1, 1}, Verbose -> False];
  data = LIESolveRegions[data, Verbose -> False];
  data = LIEExpandSeries[data, "Order" -> 4, Modulus -> modulus, "Increment" -> 2, 
    "LayerByLayer" -> True, Verbose -> False];

  data = LIEGetRelations[data, Verbose -> False, "MaxCoefDeg" -> 2];
  relations = data["Relations", "Relations"];

  Print["Matrix dimensions: ", Dimensions[relations]];
  Print[""];
  Print["(lev,deg) -> #relations:"];
  Print["         deg=0  deg=1  deg=2"];
  Do[
    levStr = "lev=" <> ToString[lev] <> "  ";
    Do[
      n = Length[relations[[lev+1, deg+1]]];
      levStr = levStr <> "  " <> ToString[n] <> "   ";
    , {deg, 0, 2}];
    Print[levStr];
  , {lev, 0, 2}];
  Print[""];

, {famName, {"bub00", "bub10", "bub11"}}];
