(* ::Package:: *)
(* VerifyExpand-SeriesVerify.wl — Series expansion substitution into IBP, truncation verification *)
(* Usage: wolframscript -file VerifyExpand-SeriesVerify.wl <famname> [order] *)
(*
    Principle:
      IBP equation (Large-Index form): Σ_α c_α(n, ν) · g[ν+α] = 0
      where g[ν+α] = J(ν+α+θn) and c_α(n, ν) are polynomials in n and ν.

      Substitute asymptotic expansion (A^(+ν) convention matching eigen-substitution):
        g(ν+α) = A^((ν+α)) · Σ_{k≥0} h^(k)(ν+α) / n^k

      where A = (A_1, ... , A_ne) are the characteristic exponents.
      NOTE: Do NOT shift ν → ν+θn in the coefficients. The operator factors
      (n+ν_a) already contain the large parameter n; shifting would double-count
      n-dependence (turning n+ν_a into 2n+ν_a). Instead, expand c_α(n, ν)
      directly in n:

        c_α(n, ν) = Σ_j c_α^(j)(ν) · n^j

      The residual at 1/n^m is then:

        Σ_{α, j}  A^(α) · c_α^(j)(ν) · h^(j+m)(ν+α)  ≡ 0  (mod p)

      where c_α^(j)(ν) = Coef_{n^j}[c_α(n, ν)],
            A^(α) = ∏ A_i^(α_i)  (α_i are integer shifts),
            h^(k)(ν+α) = shift-h^(k): ν_i → ν_i + α_i,
            h^(k) = 0 for k<0 or k>Kmax.
*)

SetDirectory[DirectoryName[$InputFileName]];
Get["VerifyExpand-Prepare.wl"];  (* shared utilities *)
Get["LIEWorkflow.wl"];           (* LIEDefineFamily *)

(* ---- Parse CLI ---- *)
famname = parseFamilyArg[];
order = If[Length[$ScriptCommandLine] >= 3,
    ToExpression[$ScriptCommandLine[[3]]],
    4
];

config = loadFamilyConfig[famname];
target = targetDir[famname];
modulus = config["Modulus"];

(* theta = (1,...,1): all nu_i → infinity at the same rate *)
theta = Table[1, {Length[config["Propagators"]]}];

Print["========================================"];
Print["VerifyExpand-SeriesVerify     fam = ", famname];
Print["  Description: ", config["Description"]];
Print["  Modulus: ", modulus];
Print["  theta = ", theta];
Print["  Order: ", order];
Print["========================================"];

(* ============================================================ *)
(* Step 1: Load expansion h^(k)                                 *)
(* ============================================================ *)
Print["\n[Step 1] Loading expansion..."];

ClearAll[$ExpansionResults];

mmaFile = findResultFile[famname, "VerifyExpansion-MMAExpansion.m"];
If[!FileExistsQ[mmaFile],
    mmaFile = findResultFile[famname, "Compare-MMAResult-" <> famname <> ".m"];
];
If[!FileExistsQ[mmaFile],
    Print["[FAIL] No expansion result found."];
    Print["  Run: wolframscript -file VerifyExpand-MMAExpand.wl ", famname];
    Exit[1];
];

raw = Get[mmaFile];
Print["  Loaded: ", mmaFile];

(* Extract expansion data + ring metadata *)
(* New format (MMAExpand): {<|"Kmax"->..., "NB"->..., "Solutions"->..., "QRingMinPoly"->..., etc. |>} *)
(* Legacy format: $ExpansionResults = {{<|...|>}} or raw list *)
expMeta = Null;
If[AssociationQ[raw[[1]]],
    expMeta = raw[[1]];
    hlist = expMeta[["Solutions"]][[1]][["H"]];
    Print["  Format: MMAExpand (NB=", expMeta["NB"], ")"];
,
    (* Handle legacy / $ExpansionResults formats *)
    If[ValueQ[$ExpansionResults] && Depth[$ExpansionResults] >= 4,
        expMeta = Quiet[$ExpansionResults[[1, 1]]];
        hlist = expMeta[["Solutions"]][[1]][["H"]];
    ];
    If[expMeta === Null && Depth[raw] >= 3 && AssociationQ[raw[[1]]],
        expMeta = raw[[1]];
        hlist = expMeta[["Solutions"]][[1]][["H"]];
    ];
    If[expMeta === Null && Depth[raw] >= 4 && ListQ[raw[[1]]] && AssociationQ[raw[[1, 1]]],
        expMeta = raw[[1, 1]];
        hlist = expMeta[["Solutions"]][[1]][["H"]];
    ];
];

If[expMeta === Null,
    Print["[FAIL] Unknown expansion format. Depth=", Depth[raw]];
    Exit[1];
];

(* Unwrap single-element lists: {{h0},{h1},...} -> {h0,h1,...} *)
If[ListQ[First[hlist]] && Length[First[hlist]] == 1 && expMeta["NB"] == 1,
    hlist = First /@ hlist
];
Print["  HList loaded. Orders: 0..", Length[hlist] - 1, "  NB=", expMeta["NB"]];

(* ============================================================ *)
(* Step 2: IBP equations                                        *)
(* ============================================================ *)
Print["\n[Step 2] IBP equations..."];

ckdata = Null;
checkpointFile = target <> "PrepareCheckpoint-" <> famname <> ".wdx";
If[FileExistsQ[checkpointFile],
    Print["  Loading from checkpoint..."];
    ckdata = Get[checkpointFile];
    If[ckdata === Null || !AssociationQ[ckdata] || !KeyExistsQ[ckdata, "Family"],
        Print["  Checkpoint invalid (", Head[ckdata], "), will re-generate..."];
        ckdata = Null;
    ];
];
If[ckdata === Null,
    Print["  Generating via LIEDefineFamily (IBP eqs only)..."];
    ckdata = LIEDefineFamily[
        config["Propagators"], config["LoopMomenta"],
        config["ExternalMomenta"], config["KinematicRules"],
        config["TopSector"],
        "Numeric" -> config["Numeric"],
        Modulus -> config["Modulus"],
        Verbose -> False
    ];
    (* Skip LIESolveRegions - SeriesVerify only needs Family.IBPEqs *)
];
ibpeqs = ckdata[["Family"]][["IBPEqs"]];
ne  = Length[ckdata[["Family"]][["AList"]]];

(* LIEFamilyDefine uses String atoms for all variables ("g", "v1", "n", etc.).
   Convert all String atoms to Symbols for pattern matching and substitution.
   Use targeted rules — blanket Replace[_, s_String :> Symbol[s]] creates
   Symbol-identity issues where Head[expr]===g is False. *)
ibpeqs = ibpeqs /. "g"[args___] :> g[args];
nuSyms = Table[Symbol["v" <> ToString[i]], {i, ne}];
nSym = Symbol["n"];
strToSymRules = Map[# -> Symbol[#] &, Join[{"n"}, Table["v" <> ToString[i], {i, ne}]]];
ibpeqs = ibpeqs /. strToSymRules;

(* Also convert String variables in h^(k) to Symbols for consistency *)
hlist = hlist /. strToSymRules;

Print["  NE = ", ne, "  #IBPEqs = ", Length[ibpeqs]];
Print["  nu-vars = ", nuSyms];

(* ============================================================ *)
(* Step 3: Quotient ring data                                    *)
(* ============================================================ *)
Print["\n[Step 3] Quotient ring setup..."];

qrNb = expMeta["NB"];
If[qrNb > 1,
    (* ---- NB > 1: general quotient ring via companion matrix ---- *)
    qrAFree = expMeta["QRingAFree"];        (* e.g. {"A"[2]} *)
    qrMinPolyList = expMeta["QRingMinPoly"]; (* e.g. {1+44856173*A[2]+A[2]^2} *)
    qrVarRule = expMeta["QRingVarRule"];     (* e.g. {"A"[1]->119616450+59808225*A[2]} *)
    qrModulus = expMeta["Modulus"];

    (* Extract the single free variable *)
    If[Length[qrAFree] =!= 1,
        Print["[FAIL] Only single free A-variable supported (VarDeg>1). Got ", qrAFree];
        Exit[1];
    ];
    freeAVar = qrAFree[[1]];  (* e.g. "A"[2] *)

    (* MinPoly: c0 + c1*A + ... + c_{d-1}*A^{d-1} + A^d = 0 *)
    minPoly = qrMinPolyList[[1]];
    deg = Exponent[minPoly, freeAVar];
    If[deg =!= qrNb,
        Print["[FAIL] MinPoly degree ", deg, " != NB=", qrNb];
        Exit[1];
    ];
    Print["  Free A: ", freeAVar, "  deg=", deg, "  mod=", qrModulus];

    (* Companion matrix M1 for multiplication by A in F_p[A]/(MinPoly) *)
    (* A^d = -c_0 - c_1*A - ... - c_{d-1}*A^{d-1} *)
    M1 = ConstantArray[0, {qrNb, qrNb}];
    Do[M1[[i + 1, i]] = 1, {i, qrNb - 1}];
    Do[M1[[i, qrNb]] = PolynomialMod[-Coefficient[minPoly, freeAVar, i - 1], qrModulus], {i, qrNb}];

    (* basisMatrices[[k+1]] = M1^k represents multiplication by A^k *)
    (* Build iteratively: M1^0 = I, M1^{k+1} = M1^k . M1 (mod p) *)
    basisMatrices = {IdentityMatrix[qrNb]};
    Do[
        AppendTo[basisMatrices, PolynomialMod[basisMatrices[[-1]] . M1, qrModulus]];
    , {qrNb - 1}];

    (* ---- General ring operations for any NB ---- *)
    (* Convert NB-vector (coeffs in monomial basis {1,A,...,A^{d-1}}) to NB*NB matrix *)
    ringElementToMatrix[vec_] := Sum[vec[[k + 1]] * basisMatrices[[k + 1]], {k, 0, qrNb - 1}];

    (* Multiply two ring elements *)
    ringMultiply[vec1_, vec2_] := PolynomialMod[ringElementToMatrix[vec1] . vec2, qrModulus];

    (* Inverse via LinearSolve: mat(x) . x^{-1} = e1 *)
    ringInverse[vec_] := LinearSolve[ringElementToMatrix[vec], UnitVector[qrNb, 1], Modulus -> qrModulus];

    (* Power with exponent in Z *)
    ringPower[vec_, exp_] := Module[{res},
        If[exp == 0, Return[UnitVector[qrNb, 1]]];
        If[exp > 0,
            res = UnitVector[qrNb, 1];
            Do[res = ringMultiply[res, vec], {exp}];
            res,
            ringPower[ringInverse[vec], -exp]
        ]
    ];

    (* Build A_i NB-vectors from VarRule *)
    freeAIdx = ToExpression[StringTake[ToString[freeAVar], {3}]];  (* "A"[2] -> 2 *)
    aVectors = Table[ConstantArray[0, qrNb], {ne}];
    aVectors[[freeAIdx]] = UnitVector[qrNb, 2];  (* A itself = {0, 1, 0, ...} *)

    Do[
        lhs = qrVarRule[[ri, 1]];  (* "A"[i] *)
        rhs = qrVarRule[[ri, 2]];  (* expression in free A *)
        iIdx = ToExpression[StringTake[ToString[lhs], {3}]];
        aVectors[[iIdx]] = Table[
            PolynomialMod[Coefficient[rhs, freeAVar, k], qrModulus],
            {k, 0, qrNb - 1}
        ];
    , {ri, Length[qrVarRule]}];

    Print["  A-vectors: ", MapIndexed[{"A[", #2[[1]], "]=", #1} &, aVectors]];

    (* aPlusShift: prod_i A_i^{shift_i} in quotient ring *)
    aPlusShift[shift_] := Module[{res = UnitVector[qrNb, 1]},
        Do[
            res = ringMultiply[res, ringPower[aVectors[[i]], shift[[i]]]];
        , {i, ne}];
        res
    ];

    (* evalH: h^(k)(nu+shift) as NB-vector of nu-polynomials *)
    evalH[k_, shift_] := Module[{rule},
        If[k < 0 || k >= Length[hlist], Return[ConstantArray[0, qrNb]]];
        rule = Thread[nuSyms -> (nuSyms + shift)];
        Map[PolynomialMod[# /. rule, qrModulus] &, hlist[[k + 1]]]
    ];

    Print["  General ring helpers defined (NB=", qrNb, ", mod=", qrModulus, ")."];
,
    (* ---- NB == 1: plain F_p (backward compatible) ---- *)
    qrModulus = modulus;

    (* Use A-values directly from the export (ARegList-aligned) *)
    (* VarRule has form: {"A"[1] -> val1, "A"[2] -> val2, ...} *)
    varRule = expMeta[["QRingVarRule"]];
    aVals = Table[
        aSymStr = StringJoin["A[", ToString[i], "]"];
        (* Look up "A"[i] -> val, convert to numeric *)
        val = Lookup[varRule, "A"[i], Null];
        If[val === Null,
            Print["[WARN] A[", i, "] not found in VarRule, assuming 0"];
            val = 0
        ];
        val
    , {i, ne}];
    Print["  A-values (F_p): ", aVals];

    (* Simple F_p helpers *)
    aPlusShiftFP[shift_] := PolynomialMod[
        Product[PowerMod[aVals[[i]], shift[[i]], modulus], {i, ne}],
        modulus
    ];

    evalHFP[k_, shift_] := Module[{rule},
        If[k < 0 || k >= Length[hlist], Return[0]];
        rule = Thread[nuSyms -> (nuSyms + shift)];
        PolynomialMod[hlist[[k + 1]] /. rule, modulus]
    ];
    Print["  F_p helpers defined (mod=", modulus, ")."];
];

(* ============================================================ *)
(* Step 4: Parse IBP equations                                  *)
(* ============================================================ *)
Print["\n[Step 4] Parsing IBP equations..."];

(*
   Each IBP eq is a sum of terms:  coef(n, nu) * g[nu + alpha]
   where g[nu+alpha] = G(nu+alpha). Extract integer shift alpha_i
   from g[v1+alpha1, v2+alpha2] by setting all nu_i -> 0.
*)
parseIBPEq[eq_] := Module[{terms, parsed, gPart, gArgs, shift, coef, grouped, result},
    terms = If[Head[eq] === Plus, List @@ eq, {eq}];
    parsed = {};
    Do[
        gPart = Cases[term, _g, {0, Infinity}];
        If[gPart === {},
            AppendTo[parsed, {Expand[term], Table[0, {ne}]}];
        ,
            gPart = gPart[[1]];
            gArgs = List @@ gPart;
            shift = Table[gArgs[[i]] /. Thread[nuSyms -> 0], {i, ne}];
            coef = Expand[term / gPart];
            AppendTo[parsed, {coef, shift}];
        ];
    , {term, terms}];

    (* Combine terms sharing the same shift *)
    grouped = GroupBy[parsed, Last];
    result = {};
    Do[
        sum = PolynomialMod[Total[First /@ grouped[shift]], modulus];
        AppendTo[result, {sum, shift}];
    , {shift, Keys[grouped]}];
    result
];

parsedEqs = parseIBPEq /@ ibpeqs;
Print["  Parsed ", Length[parsedEqs], " equations."];
If[Length[parsedEqs] > 0,
    Print["  First eq (", Length[parsedEqs[[1]]], " distinct shifts):"];
    Do[
        {c, s} = parsedEqs[[1, t]];
        Print["    shift=", s, "  coef=", Short[c]];
    , {t, Min[Length[parsedEqs[[1]]], 4]}];
];

(* ============================================================ *)
(* Step 5: Verify truncation                                    *)
(* ============================================================ *)
Print["\n[Step 5] Verifying series truncation (order=", order, ")..."];

(* Verification helpers dispatch on qrNb *)
If[qrNb > 1,
    (* ---- NB > 1: quotient ring verification ---- *)
    verifyEqAtOrder[parsedEq_, m_] := Module[{resVec, coef, shift, aqVec,
        degJ, j, cj, hEval, termVec},
        resVec = ConstantArray[0, qrNb];
        Do[
            {coef, shift} = term;
            aqVec = aPlusShift[shift];  (* NB-vector in F_p *)

            degJ = Exponent[coef, n];
            For[j = 0, j <= degJ, j++,
                cj = Coefficient[coef, n, j];
                If[cj === 0, Continue[]];
                hEval = evalH[j + m, shift];  (* NB-vector of nu-polynomials *)
                If[AllTrue[hEval, # === 0 &], Continue[]];
                (* Multiply: aqVec * cj * hEval. First cj * hEval component-wise, then ring-multiply by aqVec *)
                termVec = ringMultiply[aqVec, Map[PolynomialMod[cj * #, qrModulus] &, hEval]];
                resVec = MapThread[PolynomialMod[#1 + #2, qrModulus] &, {resVec, termVec}];
            ];
        , {term, parsedEq}];
        resVec
    ];
,
    (* ---- NB == 1: plain F_p verification ---- *)
    verifyEqAtOrder[parsedEq_, m_] := Module[{result, coef, shift, aFac,
        degJ, j, cj, hEval},
        result = 0;
        Do[
            {coef, shift} = term;
            aFac = aPlusShiftFP[shift];

            degJ = Exponent[coef, n];
            For[j = 0, j <= degJ, j++,
                cj = Coefficient[coef, n, j];
                If[cj === 0, Continue[]];
                hEval = evalHFP[j + m, shift];
                If[hEval === 0, Continue[]];
                result = result + aFac * cj * hEval;
            ];
        , {term, parsedEq}];
        PolynomialMod[Expand[result], modulus]
    ];
];

(* ---- Run verification ---- *)
totalEqs = Length[parsedEqs];
(* max n-degree in any coefficient: determines highest m we can check *)
maxDegJ = Max[Exponent[First[#], n] & /@ Flatten[parsedEqs, 1]];
maxCheck = Min[order, Length[hlist] - 1 - maxDegJ];
If[maxCheck < 0, maxCheck = 0];
Print["\n  max n-degree in coefficients: ", maxDegJ];
Print["  Verifiable range: m = 0..", maxCheck];
failedOrders = ConstantArray[0, maxCheck + 1];
failures = 0;

Do[
    eq = parsedEqs[[ei]];
    Print["\n  --- Eq ", ei, "/", totalEqs, " (", Length[eq], " shift groups) ---"];
    eqFail = False;
    Do[
        res = verifyEqAtOrder[eq, m];
        pass = If[qrNb > 1,
            AllTrue[res, # === 0 &],
            res === 0
        ];
        If[pass,
            Print["    [PASS] m=", m, " (1/n^", m, ")"],
            Print["    [FAIL] m=", m, " (1/n^", m, ")  residual = ", Short[res]];
            eqFail = True;
            failedOrders[[m + 1]] += 1;
        ];
    , {m, 0, maxCheck}];
    If[eqFail, failures++];
, {ei, totalEqs}];

(* ============================================================ *)
(* Step 6: Summary                                              *)
(* ============================================================ *)
Print["\n========================================"];
Print["SeriesVerify Summary:"];
Print["  Family: ", famname];
Print["  Order: ", order, "  theta = ", theta];
Print["  Total IBP equations: ", totalEqs];
Print["  Equations with failures: ", failures];
Do[
    If[failedOrders[[m + 1]] > 0,
        Print["    m=", m, " (1/n^", m, "): ", failedOrders[[m + 1]], " equations"]
    ];
, {m, 0, maxCheck}];

If[failures == 0,
    Print["\n[PASS] All IBP equations verified! Expansion satisfies\n",
          "       shifted IBP relations up to 1/n^", maxCheck, "."],
    Print["\n[FAIL] ", failures, " equation(s) have non-zero residual."];
];

Print["========================================"];

If[failures > 0, Exit[1]];
