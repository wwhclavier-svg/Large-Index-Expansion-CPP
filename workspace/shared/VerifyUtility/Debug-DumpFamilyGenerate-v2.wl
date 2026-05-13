(* Debug-DumpFamilyGenerate-v2.wl — Clean intermediate dumps for comparison
   Usage: wolframscript -file Debug-DumpFamilyGenerate-v2.wl bub00 *)

SetDirectory[DirectoryName[$InputFileName]];
Get["LIEUtility.wl"];
Get["LIEFamilyDefine.wl"];
Get["LIECoreAlgebra.wl"];
Get["LIERegions.wl"];
Get[ParentDirectory[DirectoryName[$InputFileName]] <> "/FamilyDatabase/FamilyDatabase.wl"];

If[Length[$ScriptCommandLine] < 2, Print["Usage: ... <famname>"]; Exit[1]];
famname = $ScriptCommandLine[[2]];
config = $FamilyDatabase[famname];

Print["=== MMA DUMP: ", famname, " ==="];

pdlist = config["Propagators"];
loopmom = config["LoopMomenta"];
extmom = config["ExternalMomenta"];
spsRep = config["KinematicRules"];
topsector = config["TopSector"];
char = config["Modulus"];

ne = Length[pdlist]; nl = Length[loopmom]; nE = Length[extmom]; nibp = nl(nl+nE);
Alist = Table["A"[i], {i, ne}];
nlist = Table["n"<>ToString[i], {i, ne}];
vlist = Table["v"<>ToString[i], {i, ne}];

Print["ne=", ne, " nl=", nl, " nE=", nE, " nibp=", nibp];
Print["char=", char];

(* === STEP 1: SP2PD === *)
Print["\n=== 1a. SP2PD ==="];
scalarproducts = Join[Union@Flatten@Outer[Times,loopmom,loopmom],Flatten@Outer[Times,loopmom,extmom]];
Print["scalarproducts=", scalarproducts];
pd2sp = Expand[pdlist]/.Thread[Rule[scalarproducts,Array["sp",ne]]]/.spsRep;
Print["pd2sp=", pd2sp];
sp2pd = LinearSolve[#[[2]],-#[[1]]+Array["z",ne]]&@CoefficientArrays[Expand[pdlist]/.Thread[Rule[scalarproducts,Array["sp",ne]]]/.spsRep,Array["sp",ne]];
Print["sp2pd=", sp2pd];
sprule = Thread[Rule[scalarproducts,sp2pd]];
Print["sprule=", sprule];

(* === STEP 2: derivL (derivative of propagators w.r.t loop momenta) === *)
Print["\n=== 1b. derivL ==="];
derivlRaw = (Expand[Expand[Outer[Times,D[#,{loopmom}],Join[loopmom,extmom]]]/.spsRep/.sprule])&/@pdlist;
Do[
  Do[
    Do[
      If[derivlRaw[[i,j,k]] =!= 0,
        Print["derivL[",i-1,"][",j-1,"][",k-1,"] = ", Expand[derivlRaw[[i,j,k]]]]
      ],
    {k, Length[derivlRaw[[i,j]]]}],
  {j, Length[derivlRaw[[i]]]}],
{i, Length[derivlRaw]}];

(* === STEP 3: genIBP (IBP identity assembly) === *)
Print["\n=== 1c. IBP equations (f-form) ==="];
mon2F[poly_,zlist_,fsymb_]/;Head[poly]=!=List:=Total[(fsymb@@Exponent[#,zlist])*(#/.(#->1&/@zlist))&/@If[Head[poly]===Plus,List@@poly,List@poly]];

genRel[j_Integer,k_Integer,ind_,derivmat_]:=Module[{eqmon,eqf,npd=Length@derivmat},
eqmon=If[j==k,"d",0]+Sum[-ind[[i]]*derivmat[[i,j,k]]/"z"[i],{i,npd}];
eqmon=Expand[eqmon/(Times@@Thread@Power["z"/@Range[npd],ind])];
eqf=mon2F[eqmon,"z"/@Range[npd],"f"]/.{"f"[vec__]:>"f"@@(-{vec})};
Return[eqf//Collect[#,_f]&];];

ibpEqsF = genRel[#[[1]],#[[2]],nlist,derivlRaw]&/@(Join@@Array[{#1,#2}&,Dimensions[derivlRaw][[2;;3]]]);
Print["#IBP equations (f-form) = ", Length[ibpEqsF]];
Do[Print["f-eq[",i-1,"] = ", ibpEqsF[[i]]], {i, Length[ibpEqsF]}];

(* === STEP 4: LargeIndexIBP (g-operator form) === *)
Print["\n=== 1d. LargeIndexIBP (g-operator) ==="];
limitSectorAll = Table[1, ne];
limitRule = Table[nlist[[i]] -> If[limitSectorAll[[i]]>0, "n"+vlist[[i]], vlist[[i]]], {i, ne}];
ibpEqsG = Take[ibpEqsF, nibp];  (* IBP only *)
ibpEqsG = ibpEqsG/.{"f"[ind__]:>"g"@@({ind}-nlist+vlist)}/.limitRule;
Print["#IBP (g-form) = ", Length[ibpEqsG]];

(* Parse each g-equation into terms: {coeff, gShift, nIdx, hasD} *)
parseGTerm[expr_] := Module[{terms,parsed},
  terms = If[Head[expr]===Plus, List@@expr, {expr}];
  parsed = Table[
    Module[{t = terms[[ti]], coeff, gshift, nIdx=0, hasD=0},
      (* Check for d-term *)
      If[!FreeQ[t, "d"],
        coeff = t/.{"d"*rest_:>rest,"d":>1} /. {c_*"g"[__]:>c, "g"[__]:>1};
        hasD = 1;
        gshift = t/.{"d"*c_*"g"[a__]:>{a}, "d"*"g"[a__]:>{a}, c_*"g"[a__]:>{a}, "g"[a__]:>{a}} /. {vlist[[j_]]:>0, "n"->0},
      (* else: n-term *)
        coeff = Coefficient[t, "n"] /. {"n"->1};
        If[coeff === 0, coeff = t/.{c_*"g"[__]:>c,"g"[__]:>1}];
        gshift = t/.{c_*"g"[a__]:>{a},"g"[a__]:>{a}} /. {vlist[[j_]]:>0, "n"->1};
        (* Find which n_i contributes *)
        Do[If[!FreeQ[t, nlist[[ii]]], nIdx = ii; Break[]], {ii, ne}];
      ];
      {coeff, gshift, nIdx, hasD}
    ], {ti, Length[terms]}];
  parsed
];

Do[
  parsed = parseGTerm[ibpEqsG[[i]]];
  Print["g-eq[", i-1, "] ", Length[parsed], " terms:"];
  Do[
    Print["  term[", ti-1, "]: coeff=", parsed[[ti,1]],
      " gShift=", parsed[[ti,2]],
      " nIdx=", parsed[[ti,3]],
      " hasD=", parsed[[ti,4]]];
  , {ti, Length[parsed]}];
, {i, Length[ibpEqsG]}];

(* === STEP 5: Subsectors === *)
Print["\n=== 1e. Subsectors ==="];
sectorlist = Normal@SparseArray[#->1&/@#,ne]&/@Subsets[Flatten@Position[topsector,1],{nl,ne}];
Print["#subsectors = ", Length[sectorlist]];
Do[Print["sub[",i-1,"] = ", sectorlist[[i]]], {i, Length[sectorlist]}];

(* === STEP 6: FTable extraction === *)
Print["\n=== 1f. FTable ==="];
(* Extract F0, F1, f2, F2D from g-operator equations *)
(* FTable[i]: F0 + Sum_j F1[j,j'] A_j + Sum_{j,k} f2[j,j',k] B[j,k] *)
(* In MMA: Coefficient[ibpeqs, g[alpha]] gives the coefficient *)
(* Here we map: g[a1,...,ane] at v=all-1 -> B[j'] where shift s_i = 2*v_i - a_i *)

(* For each equation m, extract the F-table entries *)
(* This mirrors the C++ extractFTable logic *)
Do[
  eq = ibpEqsG[[m]];
  terms = If[Head[eq]===Plus, List@@eq, {eq}];
  (* F0: d-terms, evaluated at v_i=1, n=1 *)
  f0 = 0;
  f1 = Table[0, {ne}];  (* F1[m][i] — coefficient of A_i *)
  f2d = Table[0, {ne}]; (* F2D[m][i] *)
  f2mat = Table[Table[0, {ne}], {ne}];  (* f2[m][i][j] *)

  Do[
    t = terms[[ti]];
    hasD = !FreeQ[t, "d"];
    (* Get g-shift at v=all-1 *)
    gshift = t/.{c_*"g"[a__]:>{a},"g"[a__]:>{a},"d"*c_*"g"[a__]:>{a},"d"*"g"[a__]:>{a}}/.{vlist[[j_]]:>0,"n"->1};
    If[!ListQ[gshift], Continue[]];

    If[hasD,
      (* F0 term: coeff of d*g[all-ones] at v=all-1 *)
      coeff = t/.{"d"*rest_:>rest,"d":>1}/.{c_*"g"[__]:>c,"g"[__]:>1};
      f0 += coeff,
    (* else *)
      (* n-term: extract F1, f2, F2D *)
      coeffN = Coefficient[t, "n"]/.{"n"->1};
      (* standard shift: s_i = 2*1 - gshift_i, v_i=1 *)
      s = Table[2 - gshift[[ii]], {ii, ne}];
      (* A_i has s_i < 0, B[j,i] has s_j > 0 *)
      posIdx = Select[Range[ne], s[[#]] > 0 &];
      negIdx = Select[Range[ne], s[[#]] < 0 &];
      If[Length[posIdx] == 0 && Length[negIdx] == 1,
        (* Pure A term *)
        f1[[negIdx[[1]]]] += coeffN
      ];
      If[Length[posIdx] == 1 && Length[negIdx] == 1,
        (* B[j,i] where j=posIdx, i=negIdx *)
        f2mat[[posIdx[[1]]]][[negIdx[[1]]]] += coeffN
      ];
      If[Length[posIdx] == 1 && Length[negIdx] == 0,
        (* F2D: B[j] with no A *)
        f2d[[posIdx[[1]]]] += coeffN
      ];
    ];
  , {ti, Length[terms]}];

  Print["m=", m-1, ": F0=", PolynomialMod[f0, char]];
  Do[If[f1[[i]] != 0, Print["  F1[i=", i, "]=", PolynomialMod[f1[[i]], char]]], {i, ne}];
  Do[If[f2d[[i]] != 0, Print["  F2D[i=", i, "]=", PolynomialMod[f2d[[i]], char]]], {i, ne}];
  Do[
    If[f2mat[[i]][[j]] != 0,
      Print["  f2[i=", i, "][j=", j, "]=", PolynomialMod[f2mat[[i]][[j]], char]]];
  , {i, ne}, {j, ne}];
, {m, nibp}];

(* === STEP 7: A/B equations per subsector === *)
Print["\n=== 2. A/B Equations ==="];
Do[
  sector = sectorlist[[si]];
  Print["\n--- Sector ", sector, " ---"];
  (* Apply sector limit *)
  sectorRule = Thread[vlist -> sector*"n" + vlist];
  sectorEqsG = ibpEqsG /. sectorRule;

  (* Build A/B equations — same logic as buildABEquations *)
  viVals = Table[1, {ne}];  (* v=all-1 *)

  abPolyList = {};
  Do[
    eq = sectorEqsG[[m]];
    terms = If[Head[eq]===Plus, List@@eq, {eq}];
    poly = 0;
    Do[
      t = terms[[ti]];
      hasD = !FreeQ[t, "d"];
      If[hasD, Continue[]]; (* skip d-terms *)
      (* Get coeff * n *)
      coeffN = Coefficient[t, "n"] /. {"n"->1};
      If[coeffN === 0, coeffN = t/.{c_*"g"[__]:>c,"g"[__]:>1}];
      (* Get coefficient multiplier for sector *)
      nIdx = 0;
      Do[If[!FreeQ[t, nlist[[ii]]], nIdx = ii; Break[]], {ii, ne}];
      scale = If[nIdx > 0 && sector[[nIdx]] == 1, 2, 1];
      coeff = PolynomialMod[coeffN * scale, char];

      (* Get gShift at v=all-1 *)
      gshift = t/.{c_*"g"[a__]:>{a},"g"[a__]:>{a}}/.{vlist[[j_]]:>0,"n"->1};
      If[!ListQ[gshift], Continue[]];

      (* Convert g[alpha] at v=all-1 to A/B monomial *)
      s = Table[2*viVals[[ii]] - gshift[[ii]], {ii, ne}];
      mon = 1;
      Do[
        If[s[[ii]] > 0, mon = mon*"B"[ii]^s[[ii]]];
        If[s[[ii]] < 0, mon = mon*"A"[ii]^(-s[[ii]])];
      , {ii, ne}];
      poly += coeff * mon;
    , {ti, Length[terms]}];
    poly = PolynomialMod[Expand[poly], char];
    If[poly =!= 0, AppendTo[abPolyList, poly]];
  , {m, nibp}];

  If[Length[abPolyList] > 0,
    Do[Print["AB[", i-1, "] = ", abPolyList[[i]]], {i, Length[abPolyList]}],
    Print["(no non-zero A/B equations)"]
  ];
, {si, Length[sectorlist]}];

(* === Now run the full pipeline with Verbose for region info === *)
Print["\n=== 3. Full Region Solve (LIEWorkflow) ==="];
data = LIEDefineFamily[pdlist, loopmom, extmom, spsRep, topsector,
  "Numeric" -> config["Numeric"], Modulus -> char, Verbose -> False];
data = LIESolveRegions[data, Verbose -> False];
expregdata = data["Regions"];

sectors = Sort[Keys[expregdata]];
Do[
  nRegs = Length[expregdata[sector]];
  Print["\n--- Sector ", sector, ": ", nRegs, " region(s) ---"];
  Do[
    coordRing = expregdata[[Key[sector], i, Key["CoordinateRing"]]];
    Print["Region[", i, "]:"];
    Print["  VarIndep = ", coordRing[["VarIndep"]]];
    Print["  VarDep = ", coordRing[["VarDep"]]];
    Print["  VarRule = ", coordRing[["VarRule"]]];
    Print["  MinPoly = ", coordRing[["MinPoly"]]];
    Print["  nb = ", Length[coordRing[["MonomialBasis"]]]];
    Print["  MonomialBasis = ", coordRing[["MonomialBasis"]]];

    (* FractionRule *)
    fracRule = coordRing[["FractionRule"]];
    Print["  FractionRule (", Length[fracRule], " entries):"];
    Do[Print["    ", k, " -> ", fracRule[[k]]], {k, Keys[fracRule]}];

    (* Multiplication matrices *)
    mmData = coordRing[["MultiplicationMatrixData"]];
    If[Length[mmData] > 0,
      nb = Length[coordRing[["MonomialBasis"]]];
      nMatrices = Length[mmData];
      Print["  MultiplicationMatrices: ", nMatrices, " matrices (", nb, "x", nb, ")"];
      Do[
        matName = mmData[[mi, 1]];
        matVals = mmData[[mi, 2]];
        Print["    Mat[", matName, "]:"];
        Do[Print["      ", matVals[[ri]]], {ri, Length[matVals]}];
      , {mi, nMatrices}];
    ];
  , {i, nRegs}];
, {sector, sectors}];

Print["\n=== DUMP COMPLETE ==="];
