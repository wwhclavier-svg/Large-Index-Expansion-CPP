(* test_dp323_performance.wl — Benchmark Singular vs MMA on DP323 hard sectors *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
Get["SingularCoordinateRing.wl"];
SetDirectory[$LIECPPPath <> "/verify/DP323_recompute/"];

(* Load DP323 checkpoint *)
Print["Loading DP323 checkpoint..."];
$FamilyPath = $LIECPPPath <> "/verify/DP323/";
Block[{Print = If[#[[1]] === "Warning", Print[#]&, Identity]&},
  data = Import[$FamilyPath <> "PrepareCheckpoint-DP323.wdx", "WDX"];
];
If[Head[data["Regions"]] === Join, data["Regions"] = data["Regions"][[1]]];

char = data["Config", "Modulus"];
ibpeqs = data["Family", "IBPEqs"];
Alist = data["Family", "AList"];
vlist = data["Family", "VList"];

(* Hard sectors: 37 (5 props), 38 (6 props) *)
hardSectors = {
  {0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0},  (* sector 37, 5 props *)
  {0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0}   (* sector 38, 6 props *)
};

bivarPrimeInfoSym = Symbol["LIECoreAlgebra`Private`bivarPrimeInfo"];
originalDownValues = DownValues[bivarPrimeInfoSym];

Do[
  sec = hardSectors[[i]];
  Print["\n========================================"];
  Print["Sector ", i + 36, ": ", sec, " (", Total[sec], " props)"];
  Print["========================================"];

  ibpeqslimited = LIEUtility`sectorLimitIBP[ibpeqs, sec, vlist];

  (* --- Singular version --- *)
  DownValues[bivarPrimeInfoSym] = {
    HoldPattern[bivarPrimeInfoSym[primelist_, Avar_, Bvar_, opts:OptionsPattern[]]] :>
      Module[{c = OptionValue[{opts}, Modulus], p = OptionValue[{opts}, "Parameters"],
        ls = OptionValue[{opts}, "LimitSector"], fr = OptionValue[{opts}, "FractionSymbol"]},
        SingularCoordinateRing`singularBivarPrimeInfo[primelist, Avar, Bvar,
          Modulus -> c, "Parameters" -> p, "LimitSector" -> ls, "FractionSymbol" -> fr, Verbose -> False]
      ]
  };

  timeSingular = AbsoluteTiming[
    resSingular = Check[
      LIECoreAlgebra`expRegSolve2[ibpeqslimited, Alist, vlist,
        "LimitSector" -> sec, Modulus -> char, Verbose -> False],
      $Failed
    ];
  ][[1]];
  Print["Singular: ", If[resSingular === $Failed, "FAILED", Length[resSingular] <> " region(s)"],
    ", time=", Round[timeSingular, 0.01], "s"];

  (* --- MMA version with timeout --- *)
  DownValues[bivarPrimeInfoSym] = originalDownValues;

  timeMMA = AbsoluteTiming[
    resMMA = TimeConstrained[
      LIECoreAlgebra`expRegSolve2[ibpeqslimited, Alist, vlist,
        "LimitSector" -> sec, Modulus -> char, Verbose -> False],
      120,
      $TimedOut
    ];
  ][[1]];
  Print["MMA: ", Switch[resMMA, $TimedOut, "TIMEOUT (>120s)", $Failed, "FAILED",
    _, Length[resMMA] <> " region(s), time=" <> ToString[Round[timeMMA, 0.01]] <> "s"]];

  (* Compare if both succeeded *)
  If[resMMA =!= $TimedOut && resMMA =!= $Failed && resSingular =!= $Failed,
    If[Length[resMMA] =!= Length[resSingular],
      Print["Region count mismatch!"]
      ,
      match = True;
      Do[
        If[resMMA[[j]]["VarIndep"] =!= resSingular[[j]]["VarIndep"] ||
           resMMA[[j]]["MonomialBasis"] =!= resSingular[[j]]["MonomialBasis"] ||
           resMMA[[j]]["MonomialBasisMatrix"] =!= resSingular[[j]]["MonomialBasisMatrix"],
          Print["Region ", j, " mismatch!"];
          match = False;
        ];
      , {j, Length[resMMA]}];
      If[match, Print["Results MATCH"]];
    ];
  ];
, {i, Length[hardSectors]}];

DownValues[bivarPrimeInfoSym] = originalDownValues;
Print["\nDone."];
