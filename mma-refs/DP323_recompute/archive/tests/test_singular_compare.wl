(* test_singular_compare.wl — Compare MMA vs Singular bivarPrimeInfo via expRegSolve2 *)
(* Usage: wolframscript -file test_singular_compare.wl <family> *)

$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
Get["SingularCoordinateRing.wl"];
SetDirectory[$LIECPPPath <> "/verify/DP323_recompute/"];

(* Parse family argument *)
If[Length[$ScriptCommandLine] < 2,
  famname = "SR212";
  ,
  famname = $ScriptCommandLine[[2]];
];

$FamilyPath = $LIECPPPath <> "/verify/" <> famname <> "/";
Print["Family: ", famname];

(* Load checkpoint *)
Print["Loading checkpoint..."];
cpFile = $FamilyPath <> "PrepareCheckpoint-" <> famname <> ".wdx";
If[!FileExistsQ[cpFile],
  Print["ERROR: Checkpoint not found: ", cpFile];
  Exit[1];
];
Block[{Print = If[#[[1]] === "Warning", Print[#]&, Identity]&},
  data = Import[cpFile, "WDX"];
];
If[Head[data["Regions"]] === Join, data["Regions"] = data["Regions"][[1]]];

char = data["Config", "Modulus"];
ibpeqs = data["Family", "IBPEqs"];
Alist = data["Family", "AList"];
vlist = data["Family", "VList"];
existingRegions = data["Regions"];

(* Pick a simple completed sector with 1 region *)
targetSector = SelectFirst[SortBy[Keys[existingRegions], Total],
  Length[existingRegions[#]] == 1 &,
  $Failed
];
If[targetSector === $Failed,
  Print["ERROR: No suitable 1-region sector found"];
  Exit[1];
];

Print["\n========================================"];
Print["Sector: ", targetSector, " (", Total[targetSector], " props, ", Length[existingRegions[targetSector]], " region)"];
Print["========================================"];

(* Apply sector limit *)
ibpeqslimited = LIEUtility`sectorLimitIBP[ibpeqs, targetSector, vlist];

(* --- Run 1: MMA original --- *)
Print["\n--- MMA original ---"];
timeMMA = AbsoluteTiming[
  resMMA = LIECoreAlgebra`expRegSolve2[ibpeqslimited, Alist, vlist,
    "LimitSector" -> targetSector, Modulus -> char, Verbose -> False];
][[1]];
Print["MMA result: ", Length[resMMA], " region(s), time=", timeMMA, "s"];

(* --- Run 2: Singular replacement --- *)
(* Temporarily replace bivarPrimeInfo with singularBivarPrimeInfo *)
bivarPrimeInfoSym = Symbol["LIECoreAlgebra`Private`bivarPrimeInfo"];
originalDownValues = DownValues[bivarPrimeInfoSym];

(* Define wrapper to pass options correctly *)
singularWrapper[primelist_, Avar_, Bvar_, opts:OptionsPattern[]] := Module[{
    char = OptionValue[{opts}, Modulus],
    params = OptionValue[{opts}, "Parameters"],
    limitSector = OptionValue[{opts}, "LimitSector"],
    FRsym = OptionValue[{opts}, "FractionSymbol"]
  },
  SingularCoordinateRing`singularBivarPrimeInfo[
    primelist, Avar, Bvar,
    Modulus -> char,
    "Parameters" -> params,
    "LimitSector" -> limitSector,
    "FractionSymbol" -> FRsym,
    Verbose -> False
  ]
];

DownValues[bivarPrimeInfoSym] = {
  HoldPattern[bivarPrimeInfoSym[primelist_, Avar_, Bvar_, opts:OptionsPattern[]]] :>
    singularWrapper[primelist, Avar, Bvar, opts]
};

Print["\n--- Singular replacement ---"];
timeSingular = AbsoluteTiming[
  resSingular = LIECoreAlgebra`expRegSolve2[ibpeqslimited, Alist, vlist,
    "LimitSector" -> targetSector, Modulus -> char, Verbose -> False];
][[1]];
Print["Singular result: ", Length[resSingular], " region(s), time=", timeSingular, "s"];

(* Restore original *)
DownValues[bivarPrimeInfoSym] = originalDownValues;

(* --- Compare --- *)
Print["\n========================================"];
Print["Comparison"];
Print["========================================"];

If[Length[resMMA] =!= Length[resSingular],
  Print["FAIL: Different region counts: MMA=", Length[resMMA], ", Singular=", Length[resSingular]];
  Exit[1];
];

keysToCheck = {"VarDep", "MinPoly", "VarIndep", "VarRule", "FractionRule",
  "MonomialBasis", "MonomialBasisIndex", "MonomialBasisMatrix", "VarDeg", "LimitSector"};

allPass = True;
Do[
  regMMA = resMMA[[i]];
  regS = resSingular[[i]];
  Print["\nRegion ", i, ":"];
  Do[
    key = keysToCheck[[k]];
    valMMA = regMMA[key];
    valS = regS[key];
    (* Normalize for comparison *)
    If[Head[valMMA] === List && Head[valS] === List,
      valMMA = Sort[valMMA];
      valS = Sort[valS];
    ];
    If[Head[valMMA] === Association && Head[valS] === Association,
      valMMA = KeySort[valMMA];
      valS = KeySort[valS];
    ];
    If[valMMA =!= valS,
      Print["  FAIL: ", key];
      Print["    MMA:     ", Short[valMMA, 50]];
      Print["    Singular:", Short[valS, 50]];
      allPass = False;
      ,
      Print["  PASS: ", key];
    ];
  , {k, Length[keysToCheck]}];
, {i, Length[resMMA]}];

Print["\n========================================"];
If[allPass,
  Print["ALL CHECKS PASSED"];
  Print["Speedup: ", Round[timeMMA/timeSingular, 0.1], "x"];
  ,
  Print["SOME CHECKS FAILED"];
];
Print["========================================"];
