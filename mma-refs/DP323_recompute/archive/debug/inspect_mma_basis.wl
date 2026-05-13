(* inspect_mma_basis.wl *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$BoxPath = $LIECPPPath <> "/verify/Box/";

baseline = Import[$BoxPath <> "PrepareCheckpoint-Box-MMA.wdx", "WDX"];
If[Head[baseline["Regions"]] === Join, baseline["Regions"] = baseline["Regions"][[1]]];

mmaRegions = baseline["Regions"][{1,1,1,1}];
Do[
  ring = mmaRegions[[i]]["CoordinateRing"];
  Print["MMA Region ", i, ":"];
  Print["  VarIndep=", ring["VarIndep"]];
  Print["  MonomialBasis=", ring["MonomialBasis"]];
  Print["  MonomialBasisIndex=", ring["MonomialBasisIndex"]];
  Print["  MonomialBasisMatrix keys=", Keys[ring["MonomialBasisMatrix"]]];
, {i, Length[mmaRegions]}];
