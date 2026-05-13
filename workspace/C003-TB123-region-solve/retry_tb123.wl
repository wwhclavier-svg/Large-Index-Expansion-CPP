(* retry_tb123.wl -- Retry error sectors for TB123 with Retry tier *)
$ProjectRoot = "/home/ykm/Large-Index-Expansion-CPP";
$VerifyUtilityPath = FileNameJoin[{$ProjectRoot, "verify", "VerifyUtility"}];
$FamilyDBPath = FileNameJoin[{$ProjectRoot, "verify", "FamilyDatabase", "FamilyDatabase.wl"}];
$CacheDir = FileNameJoin[{$ProjectRoot, "verify", "TB123", "cache"}];
$OutputDir = $ProjectRoot;

SetDirectory[$VerifyUtilityPath];
Print["========================================"];
Print["Retry TB123 error sectors"];
Print["========================================"];

(* Load packages *)
Get["LIEWorkflow.wl"];
Get["ScheduledRegionSolver.wl"];
Get["ExportBinary_IBPMatrix.wl"];

(* Load family *)
Get[$FamilyDBPath];
config = $FamilyDatabase["TB123"];

(* Define family *)
Print["\nDefining family..."];
data = LIEDefineFamily[
    config["Propagators"], config["LoopMomenta"],
    config["ExternalMomenta"], config["KinematicRules"],
    config["TopSector"],
    "Numeric" -> config["Numeric"],
    Modulus -> config["Modulus"],
    Verbose -> True
];
Print["  Sectors: ", Length[data["Family", "SectorList"]]];

(* Error sectors to retry *)
errorSectors = {"0111110", "1101010", "1101110", "1110110", "1111010", "1111110"};
Print["\nRetrying ", Length[errorSectors], " error sectors: ", errorSectors];

(* Package for ScheduledRegionSolve *)
familyConfig = <|
    "Family" -> data["Family"],
    "Config" -> <|"Modulus" -> config["Modulus"]|>
|>;

(* Run on error sectors only *)
Print["\nRunning ScheduledRegionSolve (error sectors only)..."];
scheduleConfig = <|
    "Tiers" -> {
        <|"Name" -> "Retry", "Timeout" -> 172800, "Order" -> "FIFO", "MaxParallel" -> 1|>
    },
    "MaxRetriesPerTier" -> 0,
    "CacheDir" -> $CacheDir,
    "OutputDir" -> $OutputDir,
    "LogDir" -> $CacheDir,
    "ScheduleMode" -> "Sequential",
    "MaxParallelJobs" -> 1,
    "CheckpointInterval" -> 1,
    "Verbose" -> True,
    "Sectors" -> errorSectors
|>;

result = ScheduledRegionSolve["TB123", familyConfig, scheduleConfig];

Print["\nRetry complete."];
status = result["Status"];
failed = result["FailedSectors"];
done = Length[Select[status, #["State"] === "done" &]];
err = Length[Select[status, #["State"] === "error" &]];
Print["  Done: ", done];
Print["  Remaining error: ", err];
If[Length[failed] > 0, Print["  Failed: ", failed]];

(* Export updated binaries *)
Print["\nExporting binary matrices..."];
If[FileExistsQ[FileNameJoin[{$OutputDir, "IBPMat_TB123.bin"}]],
    allRegions = LIERegions`ScheduledRegionResult`LoadRegionCache[$CacheDir];
    If[allRegions =!= $Failed && Length[allRegions] > 0,
        Alist = data["Family", "AList"]; ne = Length[Alist]; char = config["Modulus"];
        ExportIBPMatrixBinary`ExportBinaryIBPMatrix[
            FileNameJoin[{$OutputDir, "IBPMat_TB123.bin"}], allRegions, char];
        ExportIBPMatrixBinary`ExportBinaryRingData[
            FileNameJoin[{$OutputDir, "RingData_TB123.bin"}], allRegions, Alist, ne, char];
        Print["  IBPMat_TB123.bin updated"];
        Print["  RingData_TB123.bin updated"];
    ];
];
Print["========================================"];
