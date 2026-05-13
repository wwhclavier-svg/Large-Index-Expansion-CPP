(* debug_mbpolys.wl — Inspect raw Singular output for Box region 2 *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";
$FamilyPath = $LIECPPPath <> "/families/";

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
Get["SingularCoordinateRing.wl"];
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

(* Run only the Singular part to inspect raw output *)
ibpeqssub = sectorLimitIBP[ibpeqs, sector, vlist];
expregsub = expRegSolve2[ibpeqssub, Alist, vlist, "LimitSector" -> sector, Modulus -> char, Verbose -> False];

(* Get prime list for region 2 (the first non-trivial region) *)
primes = expregsub[[2]]["MinPoly"];  (* This won't work, need to get primes from expRegSolve2 internal *)
