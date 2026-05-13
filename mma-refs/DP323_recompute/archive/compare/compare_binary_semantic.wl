(* compare_binary_semantic.wl — Semantic comparison of binary exports *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$BoxPath = $LIECPPPath <> "/verify/Box/";

Print["=== Loading binary files ==="];

(* Load using binary import helpers from ExportBinary_IBPMatrix.wl *)
Get[$LIECPPPath <> "/verify/VerifyUtility/ExportBinary_IBPMatrix.wl"];

(* Read raw bytes for inspection *)
mmaIBP = BinaryReadList[$BoxPath <> "IBPMat_Box-MMA.bin", "Byte"];
singIBP = BinaryReadList["IBPMat_Box-Singular.bin", "Byte"];

Print["IBPMat bytes: MMA=", Length[mmaIBP], ", Sing=", Length[singIBP]];

(* Find first difference *)
diffPos = Select[Range[Min[Length[mmaIBP], Length[singIBP]]], mmaIBP[[#]] =!= singIBP[[#]] &];
Print["First 5 differences at positions: ", Take[diffPos, UpTo[5]]];

Do[
  pos = diffPos[[i]];
  Print["  pos ", pos, ": MMA=", mmaIBP[[pos]], ", Sing=", singIBP[[pos]]];
, {i, UpTo[10]}];

(* Load RingData similarly *)
mmaRing = BinaryReadList[$BoxPath <> "RingData_Box-MMA.bin", "Byte"];
singRing = BinaryReadList["RingData_Box-Singular.bin", "Byte"];

Print["\nRingData bytes: MMA=", Length[mmaRing], ", Sing=", Length[singRing]];

diffPos2 = Select[Range[Min[Length[mmaRing], Length[singRing]]], mmaRing[[#]] =!= singRing[[#]] &];
Print["First 5 differences at positions: ", Take[diffPos2, UpTo[5]]];

Do[
  pos = diffPos2[[i]];
  Print["  pos ", pos, ": MMA=", mmaRing[[pos]], ", Sing=", singRing[[pos]]];
, {i, UpTo[10]}];

Print["Done."];
