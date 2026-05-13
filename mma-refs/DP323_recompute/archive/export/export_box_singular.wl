(* export_box_singular.wl — Export Singular pipeline results to binary, compare with MMA baseline *)

$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";
$FamilyPath = $LIECPPPath <> "/families/";
$BoxPath = $LIECPPPath <> "/verify/Box/";

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
Get["SingularCoordinateRing.wl"];
Get["ExportBinary_IBPMatrix.wl"];
SetDirectory[$FamilyPath];
Get["FamilyDatabase.wl"];
SetDirectory[$LIECPPPath <> "/verify/DP323_recompute/"];

char = 179424673;
config = $FamilyDatabase["Box"];

Print["=== Running Singular pipeline on Box ==="];
workflow = LIEDefineFamily[
  config["Propagators"], config["LoopMomenta"], config["ExternalMomenta"],
  config["KinematicRules"], config["TopSector"],
  "Numeric" -> config["Numeric"], Modulus -> char
];

ibpeqs = workflow["Family", "IBPEqs"];
Alist = workflow["Family", "AList"];
vlist = workflow["Family", "VList"];
sectorList = workflow["Family", "SectorList"];

AbsoluteTiming[
  singularRes = regionsBySectors[ibpeqs, sectorList, Alist, vlist,
    Modulus -> char, Verbose -> False];
][[1]];

Print["=== Exporting Singular binary ==="];
ExportBinaryIBPMatrix["IBPMat_Box-Singular.bin", singularRes, char];
ExportBinaryRingData["RingData_Box-Singular.bin", singularRes, Alist, Length[Alist], char];

Print["=== Comparing with MMA baseline ==="];
mmaIBPMat = Import[$BoxPath <> "IBPMat_Box-MMA.bin", "Binary"];
singIBPMat = Import["IBPMat_Box-Singular.bin", "Binary"];

If[mmaIBPMat === singIBPMat,
  Print["IBPMat_Box: BYTE-IDENTICAL"],
  Print["IBPMat_Box: DIFFERENT"];
  Print["  MMA: ", Length[mmaIBPMat], " bytes"];
  Print["  Sing: ", Length[singIBPMat], " bytes"];
];

mmaRing = Import[$BoxPath <> "RingData_Box-MMA.bin", "Binary"];
singRing = Import["RingData_Box-Singular.bin", "Binary"];

If[mmaRing === singRing,
  Print["RingData_Box: BYTE-IDENTICAL"],
  Print["RingData_Box: DIFFERENT"];
  Print["  MMA: ", Length[mmaRing], " bytes"];
  Print["  Sing: ", Length[singRing], " bytes"];
];

Print["Done."];