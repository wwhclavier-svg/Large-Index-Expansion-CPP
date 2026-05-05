(* ::Package:: *)
(* VerifyRelation-SeriesVerify.wl — Internal verification of C++ LIE relations *)
(* via asymptotic expansion substitution (reverse of VerifyExpand-SeriesVerify). *)
(* Usage: wolframscript -file VerifyRelation-SeriesVerify.wl <famname> [lev] [deg] [order] *)
(*
    Principle:
      The C++ relation solver finds b_{α,β} such that:
        Σ_{α,β} b_{α,β} · (ν+θn)^β · g(ν-α; n) = 0      (∀ ν, n)

      where g(ν; n) = p(α) · Σ_{k≥0} h_k(ν) / n^k is the asymptotic expansion,
      and p(α) = ∏_i Ainv_i^{α_i} with Ainv_i from the quotient ring.

      Verification at each 1/n^m order at test ν:
        Σ_{α,β} b_{α,β} · Σ_{t=0}^{|β|} c_{β,t}(ν) · p(α) · h_{m+t}(ν-α) ≡ 0  (mod p)

      where c_{β,t}(ν) = coeff of n^{t} in (ν+θn)^β, and h_k = 0 for k > Kmax.

      This is an INTERNAL verification (no Kira needed) — it uses the same h_k
      expansion coefficients that were used to build the C++ matrix.

    STATUS (2026-05-05): PARTIALLY WORKING — PENDING FIX
      - The mathematical formula is CORRECT (verified by direct numerical reconstruction
        of the C++ equation matrix at specific ν points).
      - However, the C++ AdaptiveEquationBuilder only enforces equations up to
        |stable_order| (determined by order-stability convergence analysis).
        For bub00 (lev=1, deg=1), stable_order=0, meaning only n^1 and n^0 equations
        are constrained. Higher-order coefficients (n^{-1} and below) are NOT
        enforced by the solver, so they do NOT vanish.
      - Consequence: only m ≤ stable_order checks pass. The script currently tests
        m = 0..(kmax - max|beta|) without accounting for the stable_order bound,
        leading to false failures at higher m.
      - See verify/docs/Test-Relation.md §6 for detailed analysis.

    TODO:
      1. Read stable_order from $AllRelations entry and limit verification to m ≤ stable_order
      2. Or: extend C++ solver to enforce higher orders before exporting relations
*)

(* ---- Paths ---- *)
$scriptDir = DirectoryName[$InputFileName];
SetDirectory[$scriptDir];
$CPPVerifyRoot = ParentDirectory[$scriptDir] <> "/";

(* ---- CLI parsing ---- *)
parseArgs[] := Module[{args, fam, lev, deg, order},
    args = $ScriptCommandLine;
    If[Length[args] < 2,
        Print["Usage: wolframscript -file VerifyRelation-SeriesVerify.wl <famname> [lev] [deg] [order]"];
        Exit[1]
    ];
    fam = args[[2]];
    lev = If[Length[args] >= 3, ToExpression[args[[3]]], 2];
    deg = If[Length[args] >= 4, ToExpression[args[[4]]], 2];
    order = If[Length[args] >= 5, ToExpression[args[[5]]], 10];
    {fam, lev, deg, order}
];

{famname, lev, deg, order} = parseArgs[];

(* ---- Load family config ---- *)
dbFile = $CPPVerifyRoot <> "FamilyDatabase/FamilyDatabase.wl";
If[!FileExistsQ[dbFile], Print["ERROR: FamilyDatabase not found"]; Exit[1]];
Get[dbFile];
If[!KeyExistsQ[$FamilyDatabase, famname],
    Print["ERROR: Unknown family '", famname, "'."];
    Do[Print["  ", k], {k, Keys[$FamilyDatabase]}];
    Exit[1]
];
config = $FamilyDatabase[famname];
modulus = config["Modulus"];
ne = Length[config["TopSector"]];

targetDir = $CPPVerifyRoot <> famname <> "/";

Print["========================================"];
Print["VerifyRelation-SeriesVerify   fam = ", famname];
Print["  lev=", lev, "  deg=", deg, "  order=", order];
Print["  Modulus: ", modulus, "  NE: ", ne];
Print["========================================"];

(* ============================================================ *)
(* Step 1: Load expansion data h_k(ν)                            *)
(* ============================================================ *)
Print["\n[Step 1] Loading expansion..."];

cResultFile = targetDir <> "Compare-CPPResult-" <> famname <> ".m";
If[!FileExistsQ[cResultFile],
    Print["ERROR: Expansion file not found: ", cResultFile];
    Print["  Run: ./build/test_expandFF ", famname];
    Exit[1]
];
Get[cResultFile];
If[!ValueQ[$ExpansionResults],
    Print["ERROR: $ExpansionResults not defined in ", cResultFile];
    Exit[1]
];

expMeta = $ExpansionResults[[1, 1]];
hlistRaw = expMeta[["Solutions"]][[1]][["H"]];
nb = expMeta["NB"];
kmax = expMeta["Kmax"];
If[nb > 1,
    Print["ERROR: NB>1 quotient ring not yet supported. Got NB=", nb];
    Exit[1]
];
(* Unwrap single-element lists *)
If[ListQ[First[hlistRaw]] && Length[First[hlistRaw]] == 1,
    hlist = First /@ hlistRaw,
    hlist = hlistRaw
];
Print["  Loaded h_k, k=0..", Length[hlist]-1, "  NB=", nb];

(* ============================================================ *)
(* Step 2: Load region info for A values (theta/characteristic) *)
(* ============================================================ *)
Print["\n[Step 2] Loading region info..."];

regionFile = targetDir <> "Compare-RegionInfo-" <> famname <> ".m";
If[!FileExistsQ[regionFile],
    Print["ERROR: Region info not found: ", regionFile];
    Print["  Run: wolframscript -file Compare-FamilyGenerate.wl ", famname];
    Exit[1]
];
regionSummary = Get[regionFile];
(* regionSummary is a list of associations *)
region = regionSummary[[1]];  (* take first (top) sector *)
varRule = region[["VarRule"]];  (* {"A"[i] -> value, ...} *)
Print["  VarRule: ", varRule];

(* Extract A values using string-head pattern matching *)
aVals = Table[
    Quiet[Check[
        Replace["A"[i], varRule],
        Print["ERROR: A[", i, "] not found in VarRule"]; Exit[1]
    ]],
    {i, ne}
];
Print["  A-values = ", aVals];
(* A_inv = modular inverse *)
aInvVals = Table[PowerMod[aVals[[i]], -1, modulus], {i, ne}];
Print["  A-inv values = ", aInvVals];

(* theta = limitSector values *)
theta = region[["VarDep"]];
Print["  theta = ", theta];

(* deg_eff = max |beta|_{supp(theta)} = deg (since all theta_i != 0) *)
degEff = If[AllTrue[theta, # != 0 &], deg, 0];
Print["  degEff = ", degEff];

(* ============================================================ *)
(* Step 3: Load C++ relations (unified $AllRelations format)       *)
(* ============================================================ *)
Print["\n[Step 3] Loading C++ relations..."];

(* Search for unified relation file across known locations *)
relFileCandidates = {
    $CPPVerifyRoot <> "AllRelations_" <> famname <> "_k" <> ToString[order] <> ".m",
    $CPPVerifyRoot <> "../AllRelations_" <> famname <> "_k" <> ToString[order] <> ".m",
    targetDir <> "AllRelations_" <> famname <> "_k" <> ToString[order] <> ".m"
};
relFile = None;
Do[If[FileExistsQ[relFileCandidates[[i]]], relFile = relFileCandidates[[i]]; Break[]], {i, Length[relFileCandidates]}];
If[relFile === None,
    Print["ERROR: Relation file not found: AllRelations_", famname, "_k", ToString[order], ".m"];
    Print["  Searched: ", relFileCandidates];
    Print["  Run: ./build/test_relationFF ", famname, " ", order, " ", lev, " ", deg];
    Exit[1]
];
Get[relFile];

(* Find entry matching (lev, deg) *)
relData = None;
Do[
    If[$AllRelations[[i]]["Lev"] == lev && $AllRelations[[i]]["Deg"] == deg,
        relData = $AllRelations[[i]];
        Break[]
    ],
    {i, Length[$AllRelations]}
];
If[relData === None,
    Print["ERROR: No entry with lev=", lev, " deg=", deg, " in ", relFile];
    Print["  Available entries:"];
    Do[Print["    lev=", e["Lev"], " deg=", e["Deg"], " solutions=", e["NumSolutions"]], {e, $AllRelations}];
    Exit[1]
];
alphas = relData[["Alphas"]];   (* list of {α1, α2} vectors *)
betas = relData[["Betas"]];     (* list of {β1, β2} vectors *)
nAlpha = Length[alphas];
nBeta = Length[betas];
nSols = relData[["NumSolutions"]];
coeffs = relData[["Coefficients"]];  (* (nAlpha*nBeta) x (1+nSols) matrix *)

Print["  File: ", relFile];
Print["  Alphas: ", alphas];
Print["  Betas: ", betas];
Print["  Variables: ", nAlpha*nBeta, "  Solutions: ", nSols];

If[nSols == 0,
    Print["\n========================================"];
    Print["SeriesVerify Summary:"];
    Print["  Family: ", famname, "  lev=", lev, "  deg=", deg];
    Print["  No non-trivial relations (NumSolutions=0). Nothing to verify."];
    Print["========================================"];
    Exit[0]
];

(* ============================================================ *)
(* Step 4: Helper functions (mirror C++ buildFinalMatrix)       *)
(* ============================================================ *)

(* p(α) = ∏_i Ainv_i^{α_i} *)
pAlpha[alpha_List] := PolynomialMod[
    Product[PowerMod[aInvVals[[i]], alpha[[i]], modulus], {i, ne}],
    modulus
];

(* Evaluate h_k at ν - α *)
evalH[k_, nu_List, alpha_List] := Module[{shiftedNu, hPoly},
    If[k < 0 || k >= Length[hlist], Return[0]];
    hPoly = hlist[[k + 1]];
    shiftedNu = hPoly /. Table[
        Symbol["v" <> ToString[i]] -> nu[[i]] - alpha[[i]],
        {i, ne}
    ];
    PolynomialMod[shiftedNu, modulus]
];

(* Safe power: x^0 = 1 even when x = 0 *)
safePower[x_, 0] := 1;
safePower[x_, e_] := x^e;

(* c_{β,t}(ν) = coefficient of n^t in (ν+θn)^β *)
(* For θ_i ≠ 0: (ν_i+θ_i n)^{β_i} = Σ_{s=0}^{β_i} C(β_i,s) ν_i^s θ_i^{β_i-s} n^{β_i-s} *)
(* Contribution to n^t: s_i = β_i - t_i, Σ t_i = t *)
(* = coefficient of n^t = Σ_{Σ t_i = t} ∏_i C(β_i, β_i-t_i) ν_i^{β_i-t_i} θ_i^{t_i} *)
computeCBetaT[beta_List, nu_List, t_Integer] := Module[{total, t1, t2, c1, c2},
    If[t < 0 || t > Total[beta], Return[0]];
    total = 0;
    Do[
        t2 = t - t1;
        If[t2 < 0 || t2 > beta[[2]], Continue[]];
        c1 = Binomial[beta[[1]], beta[[1]] - t1] * safePower[nu[[1]], beta[[1]] - t1] * theta[[1]]^t1;
        c2 = Binomial[beta[[2]], beta[[2]] - t2] * safePower[nu[[2]], beta[[2]] - t2] * theta[[2]]^t2;
        total = total + c1 * c2;
    , {t1, 0, Min[beta[[1]], t]}];
    PolynomialMod[total, modulus]
];

(* Contribution of term (α, β) with coefficient b to the relation at 1/n^m order *)
termContribution[alpha_List, beta_List, bVal_, nu_List, m_Integer] := Module[
    {total, pA, t, cb, hVal},
    If[bVal == 0, Return[0]];
    total = 0;
    pA = pAlpha[alpha];
    Do[
        cb = computeCBetaT[beta, nu, t];
        If[cb == 0, Continue[]];
        hVal = evalH[m + t, nu, alpha];
        If[hVal == 0, Continue[]];
        total = total + cb * hVal;
    , {t, 0, Total[beta]}];
    PolynomialMod[pA * bVal * total, modulus]
];

(* ============================================================ *)
(* Step 5: Verify each basis relation                            *)
(* ============================================================ *)
Print["\n[Step 5] Verifying relations..."];

(* Test ν points: special + random *)
testNuPoints = {
    {0, 0}, {1, 0}, {0, 1}, {1, 1}, {2, 1}, {1, 2}, {2, 2},
    {3, 3}, {5, 5}, {10, 10}, {100, 100}
};

(* Computable m range: m from 0 to kmax - max|beta|, but also bounded *)
maxBetaSum = Max[Total /@ betas];
mMax = Min[order, kmax - maxBetaSum];
If[mMax < 0, mMax = 0];
Print["  max |beta| = ", maxBetaSum, ", verifiable m = 0..", mMax];

(* Extract basis vector k (1-indexed in Mext) *)
getBasisVector[k_Integer] := Table[
    coeffs[[a * nBeta + b + 1, k + 1]],
    {a, 0, nAlpha - 1}, {b, 0, nBeta - 1}
];

totalPass = 0;
totalFail = 0;

Do[
    basisVec = getBasisVector[sol];
    Print["\n  --- Basis ", sol, " ---"];

    basisPass = True;
    Do[
        Do[
            residual = 0;
            Do[
                Do[
                    bVal = basisVec[[a + 1, b + 1]];
                    If[bVal == 0, Continue[]];
                    residual = residual + termContribution[
                        alphas[[a + 1]], betas[[b + 1]], bVal, nu, m
                    ];
                , {b, 0, nBeta - 1}];
            , {a, 0, nAlpha - 1}];
            residual = PolynomialMod[residual, modulus];

            If[residual == 0,
                totalPass++,
                totalFail++;
                If[basisPass,
                    Print["    [FAIL] Basis ", sol, " at nu=", nu, " m=", m, " (1/n^", m, "): residual=", residual];
                    basisPass = False;
                ];
            ];
        , {m, 0, mMax}];
    , {nu, testNuPoints}];

    If[basisPass,
        Print["    [PASS] All ", Length[testNuPoints] * (mMax + 1), " checks passed"];
    ];
, {sol, 1, nSols}];

(* ============================================================ *)
(* Step 6: Summary                                               *)
(* ============================================================ *)
Print["\n========================================"];
Print["SeriesVerify Summary:"];
Print["  Family: ", famname, "  lev=", lev, "  deg=", deg];
Print["  Modulus: ", modulus];
Print["  Total checks: ", totalPass + totalFail];
Print["  Passed: ", totalPass];
Print["  Failed: ", totalFail];
Print["========================================"];

If[totalFail > 0,
    Print["\n[FAIL] ", totalFail, " verification(s) failed."];
    Exit[1],
    Print["\n[PASS] All verifications passed!"];
];
