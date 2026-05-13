(* generate_tb123_mma_baseline.wl — Generate MMA baseline for TB123 *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";
$FamilyPath = $LIECPPPath <> "/families/";
$TB123Path = $LIECPPPath <> "/verify/TB123/";

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
SetDirectory[$FamilyPath];
Get["FamilyDatabase.wl"];
SetDirectory[$LIECPPPath <> "/verify/DP323_recompute/"];

char = 179424673;
config = $FamilyDatabase["TB123"];

Print["=== Defining TB123 family ==="];
workflow = LIEDefineFamily[
  config["Propagators"], config["LoopMomenta"], config["ExternalMomenta"],
  config["KinematicRules"], config["TopSector"],
  "Numeric" -> config["Numeric"], Modulus -> char, Verbose -> True
];

ibpeqs = workflow["Family", "IBPEqs"];
Alist = workflow["Family", "AList"];
vlist = workflow["Family", "VList"];
sectorList = workflow["Family", "SectorList"];

Print["Total sectors: ", Length[sectorList]];

Print["\n=== Running MMA regionsBySectors (this may take a while) ==="];
totalTime = AbsoluteTiming[
  mmaRes = regionsBySectors[ibpeqs, sectorList, Alist, vlist,
    Modulus -> char, Verbose -> False];
][[1]];
Print["MMA total time: ", totalTime, " s"];
Print["Completed sectors: ", Length[Keys[mmaRes]]];

(* Save checkpoint *)
checkpoint = <|
  "Config" -> <|"Modulus" -> char|>,
  "Family" -> workflow["Family"],
  "Regions" -> mmaRes
|>;

Print["Saving checkpoint to ", $TB123Path <> "PrepareCheckpoint-TB123-MMA.wdx"];
Export[$TB123Path <> "PrepareCheckpoint-TB123-MMA.wdx", checkpoint, "WDX"];

Print["Done."];
