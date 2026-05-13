(* Recover-BinFiles-DP323.wl — Extract original region data from corrupted checkpoint and regenerate .bin files *)
SetDirectory[DirectoryName[$InputFileName]];
$VerifyUtilityDir = ParentDirectory[DirectoryName[$InputFileName]] <> "/VerifyUtility";

(* Block all printing during import to avoid massive log output *)
Block[{Print},
  Print["Loading checkpoint..."];
  data = Import["PrepareCheckpoint-DP323.wdx", "WDX"];
];

regions = data["Regions"];
If[Head[regions] === Join,
  Print["Checkpoint Regions is corrupted (Join). Extracting original..."];
  regions = regions[[1]];
  If[!AssociationQ[regions],
    Print["ERROR: Cannot recover - first arg is not Association"];
    Exit[1];
  ];
  Print["Recovered ", Length[Keys[regions]], " sectors, ", Total[Length/@Values[regions]], " regions"];
  ,
  Print["Regions OK: ", Length[Keys[regions]], " sectors"];
];

(* Regenerate .bin files *)
Print["Regenerating .bin files..."];
SetDirectory[$VerifyUtilityDir];
Get["ExportBinary_IBPMatrix.wl"];
SetDirectory[DirectoryName[$InputFileName]];

family = data["Family"];
Alist = family["AList"];
ne = Length[Alist];
char = data["Config", "Modulus"];

ExportIBPMatrixBinary`ExportBinaryIBPMatrix["IBPMat_DP323.bin", regions, char];
Print["  IBPMat_DP323.bin regenerated"];

ExportIBPMatrixBinary`ExportBinaryRingData["RingData_DP323.bin", regions, Alist, ne, char];
Print["  RingData_DP323.bin regenerated"];

(* Save clean checkpoint *)
data["Regions"] = regions;
Export["PrepareCheckpoint-DP323.wdx", data, "WDX"];
Print["Checkpoint cleaned and saved."];

Print["\nDone. Bin files recovered."];
