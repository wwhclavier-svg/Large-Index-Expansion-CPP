(* export_ringdata_np322.wl -- Re-export RingData for NP322 from checkpoint *)
$ProjectRoot = "/home/ykm/Large-Index-Expansion-CPP";
checkpointFile = FileNameJoin[{$ProjectRoot, "PrepareCheckpoint-NP322.wdx"}];
target = FileNameJoin[{$ProjectRoot, "verify", "NP322"}];
outputRoot = $ProjectRoot;

Print["========================================"];
Print["Export RingData for NP322"];
Print["Checkpoint: ", checkpointFile];
Print["========================================"];

(* Load checkpoint *)
Print["\n[1/3] Loading checkpoint..."];
data = Import[checkpointFile, "WDX"];
If[data === $Failed, Print["ERROR: Cannot load checkpoint"]; Exit[1]];
Print["  Checkpoint loaded."];

(* Get family data and region data *)
config = data["Config"];
char = config["Modulus"];
expregdata = data["SectorRegionData"];
If[expregdata === $Failed || Head[expregdata] =!= Association,
  Print["ERROR: No SectorRegionData in checkpoint"];
  Exit[1]
];
Print["  Sectors: ", Length[expregdata]];
Print["  Modulus: ", char];

(* Compute AList and ne from family definition *)
Alist = data["Family", "AList"];
ne = Length[Alist];
Print["  NE = ", ne];
Print["  AList length: ", Length[Alist]];

(* Export binary matrices *)
Print["\n[2/3] Exporting binary matrices..."];
Get[FileNameJoin[{$ProjectRoot, "verify", "VerifyUtility", "ExportBinary_IBPMatrix.wl"}]];

binFile = FileNameJoin[{$ProjectRoot, "IBPMat_NP322.bin"}];
ringFile = FileNameJoin[{$ProjectRoot, "RingData_NP322.bin"}];

Print["  Exporting IBPMat..."];
result1 = ExportIBPMatrixBinary`ExportBinaryIBPMatrix[binFile, expregdata, char];
If[result1 === $Failed,
  Print["  WARNING: IBPMat export failed"];
,
  Print["  IBPMat: ", binFile, " (", FileByteCount[binFile], " bytes)"];
];

Print["  Exporting RingData..."];
result2 = ExportIBPMatrixBinary`ExportBinaryRingData[ringFile, expregdata, Alist, ne, char];
If[result2 === $Failed,
  Print["  WARNING: RingData export failed"];
,
  Print["  RingData: ", ringFile, " (", FileByteCount[ringFile], " bytes)"];
];

(* Verify *)
Print["\n[3/3] Verification:"];
Print["  IBPMat_NP322.bin: ", FileByteCount[binFile], " bytes"];
Print["  RingData_NP322.bin: ", FileByteCount[ringFile], " bytes"];

If[FileByteCount[ringFile] > 100,
  Print["\n✅ RingData export successful."];
,
  Print["\n❌ RingData too small (< 100 bytes)"];
];

Print["\n========================================"];
