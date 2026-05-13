(* ::Package:: *)
(* Recompute-Missing-DP323.wl — Re-compute timeout sectors with extended timeout *)
(* Usage: wolframscript -file Recompute-Missing-DP323.wl *)
(*
    Loads the existing checkpoint, identifies sectors that were skipped
    (timeout or 0 regions), re-computes ONLY those with 1200s timeout,
    and merges results back into complete .bin files.
*)

SetDirectory[DirectoryName[$InputFileName]];
$MMAVerifyRoot = ParentDirectory[DirectoryName[$InputFileName]] <> "/";
$VerifyUtilityDir = $MMAVerifyRoot <> "VerifyUtility";

(* Temporarily switch to VerifyUtility to load packages (they use relative paths) *)
SetDirectory[$VerifyUtilityDir];
Get["LIEWorkflow.wl"];
(* Switch back to target directory for file operations *)
SetDirectory[DirectoryName[$InputFileName]];

(* ============================================================ *)
(* Paths *)
(* ============================================================ *)
target = DirectoryName[$InputFileName] <> "/";
checkpointFile = target <> "PrepareCheckpoint-DP323.wdx";
regionInfoFile = target <> "Compare-RegionInfo-DP323.m";
binFile = target <> "IBPMat_DP323.bin";
ringFile = target <> "RingData_DP323.bin";
backupDir = target <> "backup_before_recompute/";

(* ============================================================ *)
(* Step 1: Load existing checkpoint *)
(* ============================================================ *)
Print["========================================"];
Print["Recompute-Missing-DP323"];
Print["========================================"];

Print["\n[Step 1] Loading existing checkpoint: ", checkpointFile];
data = Import[checkpointFile, "WDX"];
If[!AssociationQ[data],
    Print["ERROR: Failed to load checkpoint"];
    Exit[1];
];

Print["  Status: ", data["Status"]];
family = data["Family"];
fullSectorList = family["SectorList"];
Print["  Full sector list: ", Length[fullSectorList], " sectors"];

existingRegions = data["Regions"];
completedSectors = Keys[existingRegions];
Print["  Completed sectors (with non-trivial regions): ", Length[completedSectors]];

(* ============================================================ *)
(* Step 2: Identify missing sectors *)
(* ============================================================ *)
Print["\n[Step 2] Identifying missing sectors..."];

missingSectors = Complement[fullSectorList, completedSectors];
Print["  Missing sectors: ", Length[missingSectors]];

If[Length[missingSectors] == 0,
    Print["\nAll sectors already computed! Nothing to do."];
    Exit[0];
];

(* Characterize missing sectors by prop count *)
Print["\n  Missing sector breakdown by active prop count:"];
Do[
    nProps = Total[s[[1 ;; 8]]];  (* first 8 positions are active props *)
    Print["    ", s, "  (", nProps, " active props)"];
, {s, Take[missingSectors, UpTo[20]]}];
If[Length[missingSectors] > 20,
    Print["    ... and ", Length[missingSectors] - 20, " more"];
];

(* ============================================================ *)
(* Step 3: Re-compute missing sectors with extended timeout *)
(* ============================================================ *)
Print["\n[Step 3] Re-computing missing sectors..."];
Print["  Timeout: 1200s (MMA TimeConstrained + Singular timeout)"];

ibpeqs = family["IBPEqs"];
Alist = family["AList"];
vlist = family["VList"];
char = data["Config", "Modulus"];

(* Backup existing data before modifying *)
If[!FileExistsQ[backupDir], CreateDirectory[backupDir]];
Print["\n  Backing up existing files to: ", backupDir];
CopyFile[binFile, backupDir <> "IBPMat_DP323.bin", OverwriteTarget -> True];
CopyFile[ringFile, backupDir <> "RingData_DP323.bin", OverwriteTarget -> True];
CopyFile[regionInfoFile, backupDir <> "Compare-RegionInfo-DP323.m", OverwriteTarget -> True];
Print["  Backup complete."];

(* Compute missing sectors *)
timeRegions = AbsoluteTiming[
    newRegions = LIESolveRegions[
        ibpeqs, missingSectors, Alist, vlist,
        Modulus -> char,
        "EnableFieldExtension" -> True,
        Verbose -> True,
        "Timeout" -> 1200
    ];
][[1]];

newCompleted = Keys[newRegions];
Print["\n  Newly completed sectors: ", Length[newCompleted]];
Print["  New regions: ", Total[Length /@ Values[newRegions]]];
Print["  Time: ", timeRegions, " s (", Round[timeRegions/60, 0.1], " min)"];

(* Check which are still missing *)
stillMissing = Complement[missingSectors, newCompleted];
If[Length[stillMissing] > 0,
    Print["\n  *** WARNING: ", Length[stillMissing], " sectors STILL failed:"];
    Do[Print["    ", s], {s, stillMissing}];
];

(* ============================================================ *)
(* Step 4: Merge results *)
(* ============================================================ *)
Print["\n[Step 4] Merging results..."];

(* Merge new regions into existing *)
allRegions = Join[existingRegions, newRegions];
totalSectors = Length[Keys[allRegions]];
totalRegions = Total[Length /@ Values[allRegions]];
Print["  Total sectors: ", totalSectors];
Print["  Total regions: ", totalRegions];

(* Update workflow data *)
data["Regions"] = allRegions;
data["Status"] = "RegionsSolved";

(* ============================================================ *)
(* Step 5: Export updated files *)
(* ============================================================ *)

(* 5a: Region Summary *)
Print["\n[Step 5a] Exporting updated region summary..."];
regionSummary = {};
sectors = Sort[Keys[allRegions]];
Do[
    nRegs = Length[allRegions[sector]];
    Do[
        coordRing = allRegions[[Key[sector], i, Key["CoordinateRing"]]];
        varDeg = coordRing[[Key["VarDeg"]]];
        AppendTo[regionSummary, <|
            "Sector" -> sector,
            "RegionIndex" -> i,
            "VarDeg" -> varDeg,
            "NumSolutions" -> Times @@ varDeg,
            "VarIndep" -> coordRing[[Key["VarIndep"]]],
            "VarDep" -> coordRing[[Key["VarDep"]]],
            "VarRule" -> coordRing[[Key["VarRule"]]],
            "MinPoly" -> coordRing[[Key["MinPoly"]]],
            "Nb" -> Length[coordRing[[Key["MonomialBasis"]]]]
        |>];
    , {i, nRegs}];
, {sector, sectors}];
Put[regionSummary, regionInfoFile];
Print["  Exported ", Length[regionSummary], " region(s)"];

(* 5b: Binary matrices *)
Print["\n[Step 5b] Exporting binary matrices..."];
ne = Length[Alist];
Get[$MMAVerifyRoot <> "VerifyUtility/ExportBinary_IBPMatrix.wl"];
ExportIBPMatrixBinary`ExportBinaryIBPMatrix[binFile, allRegions, char];
ExportIBPMatrixBinary`ExportBinaryRingData[ringFile, allRegions, Alist, ne, char];
Print["  ", binFile];
Print["  ", ringFile];

(* 5c: Updated checkpoint *)
Print["\n[Step 5c] Saving updated checkpoint..."];
Export[checkpointFile, data, "WDX"];
Print["  ", checkpointFile];

(* ============================================================ *)
(* Summary *)
(* ============================================================ *)
Print["\n========================================"];
Print["Recompute Complete!"];
Print["  Originally completed: ", Length[completedSectors], " sectors"];
Print["  Newly computed:       ", Length[newCompleted], " sectors"];
Print["  Total:                ", totalSectors, " sectors"];
Print["  Total regions:        ", totalRegions];
Print["  Still missing:        ", Length[stillMissing]];
Print["========================================"];
