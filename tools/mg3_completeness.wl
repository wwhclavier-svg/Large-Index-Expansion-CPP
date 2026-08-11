(* MG3c: Characteristic equation completeness check *)
(*
  For each family with AllRelations_*.m:
  - Load relations via importRelationData
  - Fix: convert "g"[...] to g[...] (importRelationData uses string head)
  - For each (lev,deg) with solutions, compute characteristic ideal
  - Record Groebner basis size + dimension
  - "Complete" = zero-dimensional ideal (dim=0)
*)

cppRoot = ExpandFileName["~/Large-Index-Expansion-CPP"];
outDir = FileNameJoin[{cppRoot, "results", "intermediate"}];
If[FileType[outDir] === None, CreateDirectory[outDir]];

families = {"bub00", "bub10", "bub11", "Tri", "Box",
            "SR212", "SR212-3m", "SR212-5m", "Penta1L"};

results = <||>;

Do[
  fam = families[[f]];
  Print["=== ", fam, " ==="];

  relFile = FileNameJoin[{cppRoot, "relations",
    "AllRelations_" <> fam <> "_k4.m"}];
  If[!FileExistsQ[relFile], Print["  No file.\n"]; Continue[]];

  relationCPP = importRelationData[cppRoot, fam, 4];
  If[Head[relationCPP] =!= List || Length[relationCPP] == 0,
    Print["  No relations.\n"]; Continue[]];

  ne = $AllRelations[[1]]["NE"];
  vlist = Array[ToExpression["v" <> ToString[#]] &, ne];
  (* Default to all-1s sector if not in data *)
  topsector = $AllRelations[[1]]["Sector"];
  If[Head[topsector] =!= List, topsector = Table[1, ne]];
  Print["  ne=", ne, ", sector=", topsector];

  famResults = <||>;
  Do[
    entry = $AllRelations[[idx]];
    If[Lookup[entry, "NumRelations", 0] == 0, Continue[]];
    lev = Lookup[entry, "Lev", 1]; deg = Lookup[entry, "Deg", 0];

    relPolys = relationCPP[[
      lev - $AllRelations[[1]]["Lev"] + 1, deg + 1]];
    relPolys = relPolys /. "g"[a__] :> g[a];  (* fix: string -> symbol head *)
    relPolys = DeleteCases[relPolys, 0];
    If[Length[relPolys] == 0, Continue[]];

    Print["  (lev=", lev, ", deg=", deg, "): ", Length[relPolys], " relations"];
    charEq = characteristicEquation[relPolys, topsector, vlist];
    Print["    GB size: ", Length[charEq]];
    Print["    GB: ", charEq];

    (* Check dimension:
       A zero-dimensional ideal's GB has each variable appearing
       with a leading term that is a pure power of that variable.
       In practice: if L[charEq, B[i]] has a univariate leading term for each i,
       the ideal is zero-dimensional. *)
    dim = If[Length[charEq] == 0, "empty",
      If[Length[charEq] < ne, "positive (>0)", 0]];
    (* More precise: count # of solutions via Singular *)
    If[Length[charEq] > 0,
      eqStr = StringRiffle[ToString[#] & /@ charEq, ","];
      varStr = StringRiffle[Table["B(" <> ToString[i] <> ")", {i, ne}], ","];
      dimStr = SingularRun[
        "ring r = 0, (" <> varStr <> "), dp;" <>
        "ideal I = " <> eqStr <> ";" <>
        "ideal G = groebner(I);" <>
        "int d = dim(G);" <>
        "return(d);"];
      dim = ToExpression[dimStr];
      Print["    dim = ", dim];
    ];

    famResults[ToString[lev] <> "," <> ToString[deg]] = <|
      "NumRelations" -> Length[relPolys],
      "GBSize" -> Length[charEq],
      "Dimension" -> dim,
      "Complete" -> (Head[dim] === Integer && dim == 0)
    |>,
    {idx, Length[$AllRelations]}
  ];

  results[fam] = <|
    "ne" -> ne, "topSector" -> topsector,
    "entries" -> famResults
  |>;
  Print[""],

  {f, Length[families]}
];

Export[FileNameJoin[{outDir, "mg3_completeness.json"}], results, "JSON"];
Print["\n-> ", FileNameJoin[{outDir, "mg3_completeness.json"}]];
Print["Done."]
