(* verify_ringdata_correctness.wl — Verify A * Ainv = I for Singular RingData *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";
$FamilyPath = $LIECPPPath <> "/families/";

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
Get["SingularCoordinateRing.wl"];
Get["ExportBinary_IBPMatrix.wl"];
SetDirectory[$FamilyPath];
Get["FamilyDatabase.wl"];
SetDirectory[$LIECPPPath <> "/verify/DP323_recompute/"];

char = 179424673;
config = $FamilyDatabase["Box"];

workflow = LIEDefineFamily[
  config["Propagators"], config["LoopMomenta"], config["ExternalMomenta"],
  config["KinematicRules"], config["TopSector"],
  "Numeric" -> config["Numeric"], Modulus -> char
];

ibpeqs = workflow["Family", "IBPEqs"];
Alist = workflow["Family", "AList"];
vlist = workflow["Family", "VList"];
sector = {1, 1, 1, 1};

singularRes = regionsBySectors[ibpeqs, {sector}, Alist, vlist,
  Modulus -> char, Verbose -> False];

ringData = ComputeRingData[singularRes[sector][[1]], Alist, Length[Alist], char];

Print["=== Verifying A * Ainv = I ==="];
Do[
  A = ringData["FlatA"][[i]];
  Ainv = ringData["FlatAinv"][[i]];
  product = PolynomialMod[A . Ainv, char];
  identity = IdentityMatrix[Length[A]];
  If[product === identity,
    Print["Matrix ", i, ": A * Ainv = I  PASS"],
    Print["Matrix ", i, ": A * Ainv = I  FAIL"];
    Print["  Product: ", product];
    Print["  Identity: ", identity];
  ];
, {i, Length[ringData["FlatA"]]}];

Print["Done."];
