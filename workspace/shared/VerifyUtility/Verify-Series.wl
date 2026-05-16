(* ::Package:: *)
(* Verify-Series.wl — Series expansion substitution verification *)
(* Usage standalone: wolframscript -file Verify-Series.wl <famname> [order] *)
(* Usage import:    Get["Verify-Series.wl"]; SeriesVerifyRun["bub00", "Order"->4] *)
(* Prerequisites: RelationMeta_<fam>.m + Compare-CPPResult-<fam>.m + ExpansionMMA_<fam>.m *)

ClearAll[SeriesVerifyRun];
SeriesVerifyRun[famname_String, opts___Rule] := Module[
    {order, projectRoot, config, modulus, allRelFiles, allRelFile, NE,
     cResultFile, expansionFile, metaFile, regimes, topRegime, aVals, aInvVals,
     seriesTotalPass, seriesTotalFail, secpos, vlistSym, varRule, details},

    order = Lookup[{opts}, "Order", 4];
    projectRoot = Lookup[{opts}, "ProjectRoot",
        Module[{base},
            base = If[$InputFileName =!= "", DirectoryName[$InputFileName], Directory[]];
            Which[
                FileExistsQ[FileNameJoin[{base, "families/FamilyDatabase.wl"}]], base,
                FileExistsQ[FileNameJoin[{base, "..", "families/FamilyDatabase.wl"}]], FileNameJoin[{base, ".."}],
                FileExistsQ[FileNameJoin[{base, "..", "..", "families/FamilyDatabase.wl"}]], FileNameJoin[{base, "..", ".."}],
                FileExistsQ[FileNameJoin[{base, "..", "..", "..", "families/FamilyDatabase.wl"}]], FileNameJoin[{base, "..", "..", ".."}],
                True, (Print["ERROR: Cannot locate project root"]; $Failed)
            ]
        ]
    ];
    If[projectRoot === $Failed, Return[$Failed]];

    (* ---- Load family config ---- *)
    If[!ValueQ[$FamilyDatabase] || !KeyExistsQ[$FamilyDatabase, famname],
        Get[FileNameJoin[{projectRoot, "families/FamilyDatabase.wl"}]]
    ];
    If[!KeyExistsQ[$FamilyDatabase, famname],
        Print["ERROR: Unknown family '", famname, "'"]; Return[$Failed]
    ];
    config = $FamilyDatabase[famname];
    modulus = config["Modulus"];

    (* ---- Load AllRelations_<fam>_k*.m ---- *)
    allRelFiles = FileNames["AllRelations_" <> famname <> "_k*.m", projectRoot];
    If[allRelFiles === {},
        allRelFiles = FileNames["AllRelations_" <> famname <> "_k*.m",
            FileNameJoin[{projectRoot, "relations"}]]
    ];
    If[allRelFiles === {},
        Print["ERROR: No AllRelations_", famname, "_k*.m found"]; Return[$Failed]
    ];
    allRelFiles = Select[allRelFiles,
        StringMatchQ[FileNameTake[#],
            "AllRelations_" <> famname <> "_k" ~~ DigitCharacter .. ~~ ".m"] &];
    extractK[fname_String] := ToExpression[
        StringCases[fname, "AllRelations_" <> famname <> "_k" ~~ k : DigitCharacter .. ~~ ".m" :> k][[1]]];
    allRelFile = Last[SortBy[allRelFiles, extractK]];
    Get[allRelFile];
    NE = $AllRelations[[1]]["NE"];

    Print["\n=== SeriesVerify: ", famname, " (order=", order, ") ==="];

    (* ---- Load expansion coefficients ---- *)
    expansionFile = FileNameJoin[{projectRoot, "verify", famname, "Compare-CPPResult-" <> famname <> ".m"}];
    If[!FileExistsQ[expansionFile],
        expansionFile = FileNameJoin[{projectRoot, "ExpansionMMA_" <> famname <> ".m"}];
    ];
    If[!FileExistsQ[expansionFile],
        Print["  SKIP: No expansion file found"]; Return[$Failed]
    ];
    Get[expansionFile];
    If[!ValueQ[$ExpansionResults],
        Print["  SKIP: $ExpansionResults not defined"]; Return[$Failed]
    ];

    (* ---- Load A_i values from RelationMeta ---- *)
    metaFile = FileNameJoin[{projectRoot, "RelationMeta_" <> famname <> ".m"}];
    If[!FileExistsQ[metaFile],
        metaFile = FileNameJoin[{projectRoot, "relations", "RelationMeta_" <> famname <> ".m"}];
    ];
    If[FileExistsQ[metaFile],
        Get[metaFile];
        If[ValueQ[$RelationMeta],
            regimes = $RelationMeta["Regimes"];
            Print["  A-values from RelationMeta (regimes=", Length[regimes], ")"];
        ,
            Print["  WARNING: $RelationMeta not defined, using hardcoded A=59808223"];
            regimes = {<|"A" -> Table[59808223, {NE}], "Ainv" -> Table[59808223, {NE}], "Sector" -> config["TopSector"]|>};
        ];
    ,
        Print["  WARNING: RelationMeta not found, using hardcoded A=59808223"];
        regimes = {<|"A" -> Table[59808223, {NE}], "Ainv" -> Table[59808223, {NE}], "Sector" -> config["TopSector"]|>};
    ];
    nReg = Length[regimes];

    (* ---- Verification ---- *)
    seriesTotalPass = 0; seriesTotalFail = 0;
    vlistSym = Table[Symbol["v" <> ToString[i]], {i, NE}];
    details = {};

    (* Convert $AllRelations (list) to Association keyed by {lev, deg} *)
    allRelAssoc = Association[{#Lev, #Deg} -> # & /@ $AllRelations];

    (* Per-regime expansion data: $ExpansionResults[[r, 1, "Solutions", 1, "H"]] *)
    nReg = Length[regimes];
    regionHList = Table[
        If[r <= Length[$ExpansionResults] && AssociationQ[$ExpansionResults[[r, 1]]],
            $ExpansionResults[[r, 1, "Solutions", 1, "H"]],
            None
        ],
        {r, nReg}
    ];
    If[nReg > 0 && regionHList[[1]] =!= None,
        Print["  Regions: ", nReg, ", KB=", Length[regimes[[1, "A"]]], ", NE=", NE];
    ];

    Do[
        Module[{key, entry, lev, deg, sols, alphas, betas, coeffMat, nCoeff,
                hlist, hExpr, relExpr0, idx, ai, bi, coeff, alpha, beta,
                betaPoly, reg, regime, aVals, aInvVals, limitSector,
                relExpr, shiftExpr, gShift, relSeries, passQ, regPass, regFail},
            key = k; {lev, deg} = key;
            entry = allRelAssoc[key];
            sols = entry["NumSolutions"];
            If[sols == 0, Continue[]];
            alphas = entry["Alphas"]; betas = entry["Betas"];
            coeffMat = entry["Coefficients"]; nCoeff = Length[coeffMat];

            (* Build relation expression: Σ coeff · n^β · j[ν-α] *)
            relExpr0 = 0; idx = 1;
            Do[
                alpha = alphas[[ai]];
                Do[
                    beta = betas[[bi]];
                    coeff = coeffMat[[idx, 2]];
                    If[coeff != 0,
                        betaPoly = Times @@ Thread[Power[vlistSym, beta]];
                        relExpr0 = relExpr0 + coeff * betaPoly * ("j"[Sequence @@ alpha])
                    ];
                    idx++;
                , {bi, Length[betas]}],
            {ai, Length[alphas]}];
            If[relExpr0 === 0, Continue[]];

            regPass = 0; regFail = 0;
            Do[
                regime = regimes[[reg]];
                aVals = regime["A"];
                aInvVals = regime["Ainv"];
                nb = regime["NB"] /. {_Missing -> 1};
                limitSector = regime["Sector"] /. {0 -> 0, 1 -> 1};

                hlist = regionHList[[reg]];
                If[hlist === None || Length[hlist] == 0, Continue[]];

                degEff = If[AnyTrue[limitSector, # =!= 0 &], deg, 0];

                If[nb == 1,
                    aVals = Flatten[aVals];
                    (* === NB=1: scalar A (existing code) === *)
                    relSeries = ConstantArray[0, order + 1];
                    Do[
                        beta = betas[[bi]]; alpha = alphas[[ai]];
                        col = (ai - 1) * Length[betas] + bi;
                        coeff = coeffMat[[col, 2]];
                        If[coeff == 0, Continue[]];
                        aPow = 1;
                        Do[
                            If[aVals[[j]] == 0,
                                If[alpha[[j]] == 0, Null, aPow *= 0],
                                aPow *= PowerMod[aVals[[j]], -alpha[[j]], modulus]
                            ],
                        {j, NE}];
                        If[aPow == 0, Continue[]];
                        betaSupp = Total[Pick[beta, limitSector, x_ /; x =!= 0]];
                        Do[
                            tVec = {t1, t2};
                            binom = Binomial[beta[[1]], t1] * Binomial[beta[[2]], t2];
                            theta1 = If[limitSector[[1]] == 0, If[beta[[1]] == t1, 1, 0],
                                If[beta[[1]] == t1, 1, limitSector[[1]]^(beta[[1]] - t1)]];
                            theta2 = If[limitSector[[2]] == 0, If[beta[[2]] == t2, 1, 0],
                                If[beta[[2]] == t2, 1, limitSector[[2]]^(beta[[2]] - t2)]];
                            thFac = theta1 * theta2;
                            If[thFac === 0, Continue[]];
                            totalT = Total[tVec];
                            Do[
                                kIdx = r + betaSupp - totalT - degEff;
                                If[kIdx < 0 || kIdx > order, Continue[]];
                                hk = hlist[[kIdx + 1]];
                                hShifted = hk /. Thread[vlistSym -> vlistSym - alpha];
                                term = coeff * aPow * binom * thFac *
                                    (vlistSym[[1]]^t1 * vlistSym[[2]]^t2) * hShifted;
                                relSeries[[r + 1]] += term;
                            , {r, 0, order}];
                        , {t1, 0, beta[[1]]}, {t2, 0, beta[[2]]}];
                    , {ai, Length[alphas]}, {bi, Length[betas]}];
                    relSeries = PolynomialMod[#, modulus] & /@ relSeries;
                ,
                    (* === NB>1: matrix A + vector h === *)
                    aMat = Table[Partition[aVals[[j]], nb], {j, NE}];
                    aInvMat = Table[Partition[aInvVals[[j]], nb], {j, NE}];
                    relSeries = ConstantArray[0, {order + 1, nb}];
                    Do[
                        beta = betas[[bi]]; alpha = alphas[[ai]];
                        col = (ai - 1) * Length[betas] + bi;
                        coeff = coeffMat[[col, 2]];
                        If[coeff == 0, Continue[]];
                        (* A^{-α} = ∏ A_i^{-α_i} using inverse matrices *)
                        aProd = IdentityMatrix[nb];
                        Do[
                            If[alpha[[j]] == 0, Continue[]];
                            aPow = PolynomialMod[
                                MatrixPower[aInvMat[[j]], alpha[[j]]], modulus];
                            aProd = PolynomialMod[aProd . aPow, modulus];
                        , {j, NE}];
                        betaSupp = Total[Pick[beta, limitSector, x_ /; x =!= 0]];
                        Do[
                            tVec = {t1, t2};
                            binom = Binomial[beta[[1]], t1] * Binomial[beta[[2]], t2];
                            theta1 = If[limitSector[[1]] == 0, If[beta[[1]] == t1, 1, 0],
                                If[beta[[1]] == t1, 1, limitSector[[1]]^(beta[[1]] - t1)]];
                            theta2 = If[limitSector[[2]] == 0, If[beta[[2]] == t2, 1, 0],
                                If[beta[[2]] == t2, 1, limitSector[[2]]^(beta[[2]] - t2)]];
                            thFac = theta1 * theta2;
                            If[thFac === 0, Continue[]];
                            totalT = Total[tVec];
                            Do[
                                kIdx = r + betaSupp - totalT - degEff;
                                If[kIdx < 0 || kIdx > order, Continue[]];
                                hk = hlist[[kIdx + 1]];
                                (* For NB>1, hk is nb-vector; substitute ν → ν-α *)
                                hVec = hk /. Thread[vlistSym -> vlistSym - alpha];
                                (* aProd · hVec gives nb-vector *)
                                aTimesH = aProd . hVec;
                                Do[
                                    term = coeff * binom * thFac *
                                        (vlistSym[[1]]^t1 * vlistSym[[2]]^t2) * aTimesH[[j]];
                                    relSeries[[r + 1, j]] += term;
                                , {j, nb}];
                            , {r, 0, order}];
                        , {t1, 0, beta[[1]]}, {t2, 0, beta[[2]]}];
                    , {ai, Length[alphas]}, {bi, Length[betas]}];
                    relSeries = PolynomialMod[#, modulus] & /@ Flatten[relSeries];
                ];
                passQ = AllTrue[relSeries, # === 0 &];
                If[passQ, regPass++; seriesTotalPass++, regFail++; seriesTotalFail++];
                Print["  (lev=", lev, ", deg=", deg, ", reg=", reg,
                    ") Series: ", If[passQ, "PASS", "FAIL"]];
            , {reg, nReg}];

            AppendTo[details, <|"Lev" -> lev, "Deg" -> deg,
                "Regions" -> nReg, "Pass" -> regPass, "Fail" -> regFail|>];
        ]
    , {k, Keys[allRelAssoc]}];

    Print["\n  SeriesVerify: ",
        If[seriesTotalPass + seriesTotalFail > 0,
            ToString[seriesTotalPass] <> "/" <> ToString[seriesTotalPass + seriesTotalFail] <>
            If[seriesTotalFail > 0,
                " (" <> ToString[seriesTotalFail] <> " failed)", " PASS"],
            "N/A"]];

    <|"Pass" -> seriesTotalPass, "Fail" -> seriesTotalFail,
      "Total" -> seriesTotalPass + seriesTotalFail,
      "Details" -> details|>
];

(* ---- Standalone runner: only when called as a script, not imported ---- *)
If[!TrueQ[$LoadedByBladeFamily] && StringMatchQ[$ScriptCommandLine[[1]], ___ ~~ "Verify-Series.wl"],
    If[Length[$ScriptCommandLine] >= 2,
        famname = $ScriptCommandLine[[2]],
        famname = "bub00"
    ];
    order = If[Length[$ScriptCommandLine] >= 3,
        ToExpression[$ScriptCommandLine[[3]]],
        4
    ];
    result = SeriesVerifyRun[famname, "Order" -> order];
    If[result === $Failed, Exit[1]];
    Exit[If[result["Fail"] > 0, 1, 0]];
];
