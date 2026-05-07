(* VerifyAllRelations.wl — v3: Clean finite-field verification *)
p = 179424673;

Print["Modulus p = ", p];
Print[""];

(* Step 1: Load AllRelations *)
table = Get["AllRelations_bub00_k14_Table.m"];
Print["AllRelations loaded: ", Length[Keys[table]], " config(s)"];
Print[""];

(* Step 2: Build IBP reduction in F_p at GenericPoint s=3, d=1/3 *)
(* Compute modular inverses *)
inv3 = PowerMod[3, -1, p];     (* 1/3 mod p *)
inv9 = PowerMod[9, -1, p];     (* 1/9 mod p *)
invS = PowerMod[3, -1, p];     (* 1/s mod p, s=3 *)

dMod = inv3;  (* d = 1/3 *)
sMod = 3;

Print["1/3 mod p = ", inv3];
Print["1/9 mod p = ", inv9];
Print[""];

(* Reduction coefficients: G[a,b] = r[a][b] * G[1,1] (mod p) *)
r = ConstantArray[0, {4, 4}];  (* 1-indexed, up to 3 *)

(* Use exact rational arithmetic then convert to F_p *)
toMod[rat_Rational] := Mod[Numerator[rat] * PowerMod[Denominator[rat], -1, p], p];
toMod[int_Integer] := Mod[int, p];

(* NeatIBP IBP relations evaluated at s=3, d=1/3 *)
(* Rel 1: (3-d)G[1,1] + s G[2,1] = 0  -> G[2,1] = -(3-d)/s G[1,1] *)
r[[2,1]] = toMod[-(3 - 1/3)/3];  (* -(8/3)/3 = -8/9 *)

(* Rel 2: (-6+2d)G[1,1] - 2s G[1,2] = 0 -> G[1,2] = -(-6+2d)/(2s) G[1,1] *)
r[[1,2]] = toMod[(-6 + 2/3)/(2*3)];  (* -8/9 *)

(* Rel 3: (4-d)G[2,1] + 2s G[3,1] = 0 -> G[3,1] = -(4-d)/(2s) G[2,1] *)
r[[3,1]] = toMod[-(4 - 1/3)/(2*3)] * r[[2,1]];

(* Rel 4: (-8+2d)G[1,2] - 4s G[1,3] = 0 -> G[1,3] = (-8+2d)/(4s) G[1,2] *)
r[[1,3]] = toMod[(-8 + 2/3)/(4*3)] * r[[1,2]];

(* Better: do all arithmetic symbolically then convert *)
r2 = ConstantArray[0, {4, 4}];
r2[[2,1]] = -(3 - 1/3)/3 // Together;       (* = -8/9 *)
r2[[1,2]] = (-6 + 2/3)/(2*3) // Together;   (* = -8/9 *)
r2[[3,1]] = -(4 - 1/3)/(2*3) * r2[[2,1]] // Together;   (* = 44/81 *)
r2[[1,3]] = (-8 + 2/3)/(4*3) * r2[[1,2]] // Together;   (* = -44/81 *)
r2[[2,2]] = -(-12 + 2/3)/(4 - 1/3) * r2[[3,1]] // Together;  (* = 136/81 *)
r2[[2,3]] = -((7 - 1/3)*r2[[1,3]] + r2[[2,2]])/3 // Together;
r2[[3,2]] = r2[[2,3]];
r2[[3,3]] = -((8 - 1/3)*r2[[2,3]] + 2*r2[[3,2]])/(2*3) // Together;

rMod = ConstantArray[0, {4, 4}];
Do[rMod[[a,b]] = toMod[r2[[a,b]]], {a,1,3}, {b,1,3}];
rMod[[1,1]] = 1;

Print["Reduction coefficients:"];
Print["  G[a,b] = coeff * G[1,1] (mod p)"];
Do[
    sym = If[rMod[[a,b]] > p/2, rMod[[a,b]] - p, rMod[[a,b]]];
    Print["  G[", a, ",", b, "] = ", rMod[[a,b]], "  (symmetric: ", sym, ")"];
, {a, 1, 3}, {b, 1, 3}];
Print[""];

(* Quick cross-check using exact rational arithmetic *)
Print["Cross-check (exact rational):"];
Do[
    Print["  G[", a, ",", b, "] = ", r2[[a,b]], " * G[1,1]"];
, {a, 1, 3}, {b, 1, 3}];
Print[""];

(* Verify the IBP relations are satisfied *)
Print["Verifying IBP relations in F_p:"];
ibp1 = toMod[3 - 1/3] + Mod[sMod * rMod[[2,1]], p];
ibp2 = toMod[-6 + 2/3] - Mod[2*sMod * rMod[[1,2]], p];
Print["  IBP1: ", toMod[3 - 1/3], " + ", sMod, "*", rMod[[2,1]], " = ", ibp1];
Print["  IBP2: ", toMod[-6 + 2/3], " - 2*", sMod, "*", rMod[[1,2]], " = ", ibp2];
Print[""];

(* ===================================================== *)
(* Step 3: Verification *)
(* ===================================================== *)

parseRelStr[relStr_String] := Module[{symbStr},
    symbStr = StringReplace[relStr, {"j" -> "G", " = 0" -> ""}];
    ToExpression[symbStr]
];

checkAtPoint[expr_, nu1_, nu2_] := Module[
    {e, coeffAccum, a, b, coeffVal, term, total},
    e = expr /. {nu1 -> nu1Val, nu2 -> nu2Val};

    (* The expression is of form: Sum coeff * G[a,b]
       Replace each G[a,b] with its reduction coefficient * G[1,1] *)
    (* We just accumulate the coefficient of G[1,1] *)
    coeffAccum = 0;

    (* Expand the expression *)
    e = Expand[e];

    (* If it's a sum, extract terms *)
    terms = If[Head[e] === Plus, List @@ e, {e}];

    Do[
        term = terms[[i]];
        (* Get the G[a,b] part *)
        gpart = Cases[term, G[a_, b_] :> {a, b}, Infinity];
        If[Length[gpart] == 0,
            (* Constant term - should not happen in a homogeneous relation *)
            If[term =!= 0, Return[{False, "Unexpected constant term: " <> ToString[term]}]];
            Continue[]
        ];

        {a, b} = gpart[[1]];

        (* Coefficient in front of G[a,b] *)
        If[a <= 0 || b <= 0,
            (* Massless tadpole vanishes - skip *)
            Continue[]
        ];

        coeffVal = term /. G[_, _] -> 1;

        (* If a>3 or b>3, we can't reduce it *)
        If[a > 3 || b > 3,
            Return[{False, "No IBP reduction for G[" <> ToString[a] <> "," <> ToString[b] <> "]"}]
        ];

        coeffAccum = Mod[coeffAccum + Mod[coeffVal * rMod[[a,b]], p], p];
    , {i, Length[terms]}];

    If[coeffAccum == 0,
        Return[{True, ""}],
        Return[{False, "G[1,1] coefficient = " <> ToString[coeffAccum] <> " (mod p)"}]
    ]
];

(* Sub in nu values *)
checkRelation[relStr_, nu1Val_, nu2Val_] := Module[
    {expr, e},
    expr = parseRelStr[relStr];
    e = expr /. {nu1 -> nu1Val, nu2 -> nu2Val};
    (* Expand the expression so we can extract coefficients *)
    e = Expand[e];

    coeffAccum = 0;
    terms = If[Head[e] === Plus, List @@ e, {e}];

    Do[
        term = terms[[i]];
        gpart = Cases[term, G[a_, b_] :> {a, b}, Infinity];
        If[Length[gpart] == 0, Continue[]];
        {a, b} = gpart[[1]];
        If[a <= 0 || b <= 0, Continue[]];
        If[a > 3 || b > 3, Return[{False, "G[" <> ToString[a] <> "," <> ToString[b] <> "] unreducible"}]];
        coeffVal = term /. G[_, _] -> 1;
        coeffAccum = Mod[coeffAccum + Mod[coeffVal * rMod[[a,b]], p], p];
    , {i, Length[terms]}];

    If[coeffAccum == 0, {True, ""}, {False, "G[1,1] coeff = " <> ToString[coeffAccum] <> " (mod p)"}]
];

testPoints = {{2,2}, {3,2}, {2,3}, {3,3}};

Print["========================================"];
Print["Verification"];
Print["========================================\n"];

nTested = 0; nPassed = 0;

Do[
    key = Keys[table][[k]];
    {lev, deg} = key;
    relList = table[key];
    Print["Lev=", lev, ", Deg=", deg, " (", Length[relList], " relations):"];

    Do[
        relStr = relList[[r]];
        nTested++;
        allPassed = True;

        Do[
            {nu1, nu2} = pt;
            {ok, msg} = checkRelation[relStr, nu1, nu2];
            If[!ok, Print["  FAIL (", nu1, ",", nu2, "): ", msg]; allPassed = False];
        , {pt, testPoints}];

        If[allPassed, nPassed++; Print["  Rel[", r, "]: PASS"]];
    , {r, Length[relList]}];
    Print[""];
, {k, Length[Keys[table]]}];

Print["========================================"];
Print["Summary: ", nPassed, "/", nTested, " passed"];
Print["========================================"];
