(* export_bub00_singular.wl — Generate bub00 binary via Singular pipeline *)
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
config = $FamilyDatabase["bub00"];

workflow = LIEDefineFamily[
  config["Propagators"], config["LoopMomenta"], config["ExternalMomenta"],
  config["KinematicRules"], config["TopSector"],
  "Numeric" -> config["Numeric"], Modulus -> char
];

ibpeqs = workflow["Family", "IBPEqs"];
Alist = workflow["Family", "AList"];
vlist = workflow["Family", "VList"];
sectorList = workflow["Family", "SectorList"];

Print["Running Singular pipeline for bub00..."];
time = AbsoluteTiming[
  singularRes = regionsBySectors[ibpeqs, sectorList, Alist, vlist,
    Modulus -> char, Verbose -> False];
][[1]];
Print["Singular time: ", time, " s"];

(* Export binary files with -Singular suffix *)
ExportBinaryIBPMatrix["IBPMat_bub00-Singular.bin", singularRes, char];
ExportBinaryRingData["RingData_bub00-Singular.bin", singularRes, Alist, Length[Alist], char];

Print["Done."];
