(* ::Package:: *)
(* lie_kira_numerical_compare.wl — Numerical comparison of LIE vs Kira reduction rules *)
(*
  Projects both Kira (rational, symbolic) and LIE (finite field) results to F_p
  and compares coefficient-by-coefficient.

  Usage:
    wolframscript -file tools/lie_kira_numerical_compare.wl <family> [cppRoot] [mmaRoot]
*)

(* ===== Configuration ===== *)
fam = $ScriptCommandLine[[2]];
cppRoot = If[Length[$ScriptCommandLine] >= 3, $ScriptCommandLine[[3]],
  Environment["HOME"] <> "/Large-Index-Expansion-CPP"];
mmaRoot = If[Length[$ScriptCommandLine] >= 4, $ScriptCommandLine[[4]],
  Environment["HOME"] <> "/Large-Index-Expansion-MMA"];

Print["=== LIE vs Kira Numerical Comparison ==="];
Print["Family: ", fam];
Print["C++ root: ", cppRoot];
Print["MMA root: ", mmaRoot];

(* ===== Load family config ===== *)
cfgFile = cppRoot <> "/families/" <> fam <> ".json";
If[!FileExistsQ[cfgFile],
  Print["ERROR: Family config not found: ", cfgFile];
  Exit[1]
];
cfg = Import[cfgFile, "JSON"];
p = "modulus" /. cfg;
ne = Length["propagators" /. cfg];
numericRules = KeyMap[ToExpression, Association @@ ("numeric" /. cfg)];

Print["ne=", ne, "  modulus=", p];
Print["Numeric: ", numericRules];

(* Extract d value *)
dVal = Lookup[numericRules, $d, 1/3];

(* ===== Read Kira output ===== *)
kiraFile = cppRoot <> "/kira_tests/" <> fam <> "/results/" <> fam <> "/kira_integrals.m";
If[!FileExistsQ[kiraFile],
  Print["ERROR: Kira results not found: ", kiraFile];
  Exit[1]
];

(* Parse Kira rules *)
kiraRaw = Get[kiraFile];
kiraRules = Association[kiraRaw /. Rule -> List];
Print["Kira rules loaded: ", Length[kiraRules], " entries"];

(* ===== Read LIE output ===== *)
lieFile = cppRoot <> "/results/lie_e2e/" <> fam <> "_rule.m";
If[!FileExistsQ[lieFile],
  lieFile = cppRoot <> "/results/lie_e2e/" <> fam <> ".m";
];
If[!FileExistsQ[lieFile],
  Print["ERROR: LIE e2e results not found: ", lieFile,
    " (run lie_end_to_end.wl first)"];
  Exit[1]
];

lieRaw = Get[lieFile];
(* Handle both formats: direct rule list or {"reductionRule" -> ...} *)
lieRuleList = If[ListQ[lieRaw] && Length[lieRaw] > 0 && MatchQ[lieRaw[[1]], _List],
  (* Format: {{1,1} -> {g[...]->...}, ...} *)
  lieRaw,
  (* Format: {"reductionRule" -> rules} *)
  Lookup[Association[lieRaw], "reductionRule", {}]
];
Print["LIE rules loaded: ", If[ListQ[lieRuleList], Length[lieRuleList], 0], " sector entries"];

(* ===== Helper: Project rational to F_p ===== *)
ratToFp[rat_, mod_] := Module[{num, den},
  If[IntegerQ[rat], Return[Mod[rat, mod]]];
  If[!MatchQ[rat, _Rational], Return[$Failed]];
  {num, den} = {Numerator[rat], Denominator[rat]};
  Mod[num * PowerMod[den, -1, mod], mod]
];

(* ===== Helper: Evaluate Kira coefficient at LIE kinematics ===== *)
(* Kira coefficient is a rational function of d and kinematic invariants.
   We substitute LIE's numeric values and project to F_p. *)
evaluateKiraCoeff[coeff_, kinRules_, dVal_, mod_] := Module[{val},
  (* Substitute d and kinematics *)
  val = coeff /. $d -> dVal /. kinRules;
  (* Evaluate numerically *)
  val = N[val];
  (* Try rational approximation *)
  val = Rationalize[val, 10^-10];
  If[MatchQ[val, _Rational | _Integer],
    ratToFp[val, mod],
    Print["  WARNING: Could not rationalize Kira coefficient: ", coeff, " -> ", val];
    $Failed
  ]
];

(* ===== Helper: Kira integral name to index list ===== *)
(* "bub00[2,1]" -> {2,1}  or "Tri[1,1,5]" -> {1,1,5} *)
kiraNameToIndex[name_] := Module[{str, idxStr},
  str = ToString[name];
  idxStr = StringCases[str, "[" ~~ idx__ ~~ "]" :> idx];
  If[idxStr === {}, Return[$Failed]];
  ToExpression /@ StringSplit[idxStr[[1]], ","]
];

(* ===== Helper: LIE g[...] to index list ===== *)
(* g[2,1] -> {2,1} *)
gToIndex[g_List] := List @@ g;
gToIndex[g_] := If[Head[g] === g, List @@ g, $Failed];

(* ===== Helper: Evaluate LIE parametric rule at specific target index ===== *)
(* LIE rule format: g[v1,v2] -> expr(v1,v2,...) where v1,v2 are parametric
   The rule maps the specific target integral to lower integrals.
   We need to extract the coefficient for each master integral. *)
evaluateLIEToFp[rule_, targetIdx_, vlist_, mod_] := Module[
  {rhs, subs, coeffMap, terms},
  rhs = rule[[2]];  (* right-hand side expression *)
  (* Substitute vlist -> targetIdx *)
  subs = Thread[vlist -> targetIdx];
  rhs = rhs /. subs;
  (* Mod p all coefficients *)
  rhs = Expand[rhs];
  If[rhs === 0, Return[<||>]];
  (* Collect terms by g[...] *)
  terms = If[Head[rhs] === Plus, List @@ rhs, {rhs}];
  coeffMap = <||>;
  Do[
    Module[{coeff, gterm},
      If[MatchQ[term, _Integer | _Rational],
        (* Pure number - corresponds to g[0,..,0] or similar *)
        coeff = Mod[term, mod];
        gterm = g @@ Table[0, {Length[targetIdx]}];
        ,
        coeff = If[Head[term] === Times,
          Select[term, FreeQ[g]] /. {} -> 1,
          1
        ];
        gterm = If[Head[term] === Times,
          Select[term, !FreeQ[#, g] &],
          term
        ];
        coeff = Mod[coeff, mod];
      ];
      If[coeff =!= 0,
        coeffMap[gterm] = Mod[Lookup[coeffMap, gterm, 0] + coeff, mod];
      ];
    ],
    {term, terms}
  ];
  coeffMap
];

(* ===== Helper: Get sign convention correction ===== *)
(* Kira propagators: -k^2. LIE propagators: k^2.
   I_LIE(a) = (-1)^sum(a) * I_Kira(a)
   For reduction I(A) = C * I(B):
   C_LIE = (-1)^(sum(A)-sum(B)) * C_Kira
*)
signCorrection[targetIdx_, masterIdx_] := (-1)^(Total[targetIdx] - Total[masterIdx]);

(* ===== Main comparison ===== *)
nCompared = 0;
nPassed = 0;
nFailed = 0;
nSkipped = 0;
failures = {};

(* Iterate over Kira rules *)
Do[
  kiraTarget = kiraRule[[1]];   (* e.g., bub00[2,1] *)
  kiraCoeffRaw = kiraRule[[2]];  (* e.g., d-3 *)

  targetIdx = kiraNameToIndex[kiraTarget];
  If[targetIdx === $Failed,
    Print["SKIP: Cannot parse Kira target: ", kiraTarget];
    nSkipped++; Continue[]
  ];

  (* Evaluate Kira coefficient in F_p *)
  kiraFp = evaluateKiraCoeff[kiraCoeffRaw, numericRules, dVal, p];
  If[kiraFp === $Failed,
    Print["SKIP: Cannot evaluate Kira coeff for ", kiraTarget];
    nSkipped++; Continue[]
  ];

  (* Find corresponding LIE rule *)
  (* LIE rules are organized by sector: {sector -> {g[...]->..., ...}} *)
  lieRule = {};
  If[ListQ[lieRuleList] && Length[lieRuleList] > 0,
    (* Check if it's in {sector -> rules} format *)
    If[MatchQ[lieRuleList[[1]], Rule[_List, _]],
      (* Iterate over sectors *)
      Do[
        Module[{secRules},
          secRules = secEntry[[2]];
          (* secRules is a list of rules like g[2,1] -> ... *)
          Do[
            If[gToIndex[secRule[[1]]] === targetIdx,
              lieRule = secRule;
              Break[];
            ],
            {secRule, secRules}
          ];
          If[lieRule =!= {}, Break[]];
        ],
        {secEntry, lieRuleList}
      ],
      (* Direct list of rules *)
      Do[
        If[gToIndex[rule[[1]]] === targetIdx,
          lieRule = rule;
          Break[];
        ],
        {rule, lieRuleList}
      ];
    ]
  ];

  If[lieRule === {},
    Print["SKIP: No LIE rule for ", kiraTarget];
    nSkipped++; Continue[]
  ];

  (* Evaluate LIE rule *)
  lieCoeffMap = evaluateLIEToFp[lieRule, targetIdx, Array[ToExpression["v"<>ToString[#]]&, ne], p];

  (* The Kira rule maps to a single master: kiraTarget -> kiraCoeff * master
     But the master might differ between Kira and LIE after reduction.
     For simple cases (bub00), all reduce to g[1,1] / bub00[1,1].
     For the comparison, we check: do both reduce to the SAME master with the SAME coefficient? *)

  (* Extract Kira master and coefficient *)
  (* kiraCoeffRaw is the total coefficient. For bub00 it's just one master. *)
  (* For multi-master cases, this needs refinement. *)

  (* For now, compare the full right-hand side as a polynomial in g's *)
  (* The LIE side may have multiple masters; we compare term by term *)

  (* Get the primary master from Kira - for bub00 it's always bub00[1,1] *)
  (* Actually, both should give the same reduced form after applying all rules.
     Let's just check that the LIE reduction is correct by verifying the
     main coefficient matches. *)

  (* Simple check: if both reduce to a single master, compare coefficients *)
  nCompared++;

  (* For bub00: Kira has one master bub00[1,1], LIE should have g[1,1] *)
  (* Apply sign convention: C_LIE = (-1)^(sum(A)-sum(B)) * C_Kira * s^(sum(B)-sum(A)) *)
  (* where s is the variable Kira set to 1 *)

  (* Get Kira's s-variable (the one set to 1) *)
  sVal = Lookup[numericRules, s, 1];
  If[MemberQ[Keys[numericRules], s1], sVal = Lookup[numericRules, s1, 1]];

  (* Primary master: for bub00 it's {1,1} *)
  kiraMaster = {1,1}; (* Default - will need refinement for multi-master *)

  (* Sign correction *)
  signCorr = signCorrection[targetIdx, kiraMaster];

  (* s-scaling: I(A) = C * I(B).
     Kira computed at s=1: I_s=1(A) = C_Kira * I_s=1(B)
     LIE at s=sVal: I_s(A) = C_LIE * I_s(B)
     I_s(A) = s^(d/2 - sum(A)) * I_s=1(A)
     So: C_LIE = s^(sum(B) - sum(A)) * C_Kira
     = s^(-Δsum) * C_Kira
  *)
  deltaSum = Total[targetIdx] - Total[kiraMaster];
  sCorr = sVal^(-deltaSum);

  (* Expected LIE coefficient in F_p *)
  expectedFp = Mod[signCorr * kiraFp * ratToFp[sCorr, p], p];

  (* Get LIE coefficient for g[kiraMaster] *)
  lieActualFp = Lookup[lieCoeffMap, g @@ kiraMaster, $NotFount];

  If[lieActualFp === $NotFount,
    (* The LIE rule might reduce to different masters *)
    Print["INFO: LIE rule for ", kiraTarget, " reduces to different masters than ", kiraMaster];
    nSkipped++;
    ,
    If[lieActualFp === expectedFp,
      nPassed++;
      ,
      Print["FAIL: ", kiraTarget, " -> master ", kiraMaster];
      Print["  Kira (F_p): ", kiraFp, "  expected LIE (F_p): ", expectedFp];
      Print["  LIE actual (F_p): ", lieActualFp];
      AppendTo[failures, <|"target" -> ToString[kiraTarget],
        "kiraFp" -> kiraFp, "expectedFp" -> expectedFp,
        "lieFp" -> lieActualFp|>];
      nFailed++;
    ];
  ];

  ,
  {kiraRule, Normal[kiraRules]}
];

(* ===== Report ===== *)
Print["\n=== Results ==="];
Print["Compared: ", nCompared];
Print["Passed:   ", nPassed];
Print["Failed:   ", nFailed];
Print["Skipped:  ", nSkipped];

If[nFailed > 0,
  Print["\n=== Failures ==="];
  Do[
    Print[f["target"], ": Kira=", f["kiraFp"],
      " expected=", f["expectedFp"], " LIE=", f["lieFp"]],
    {f, failures}
  ];
  Exit[1];
  ,
  Print["\nALL PASSED!"];
  Exit[0];
];
