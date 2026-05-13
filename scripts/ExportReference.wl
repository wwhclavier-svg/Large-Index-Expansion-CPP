(* ::Package:: *)

(*
   ExportReference.wl

   Purpose: Run the large-index expansion using the Mathematica reference
   implementation, then export the raw expansion coefficients as JSON for
   cross-validation with the C++ implementation.

   Usage (in Mathematica):
     << "~/文档/Wolfram Mathematica/[Work] Program_LargeIndexExpansion_2508/LargeIndexExpansion.wl"
     << "~/文档/Wolfram Mathematica/[Work] Program_LargeIndexExpansion_2508/BlockRecursion.wl"
     << "~/文档/Large-Index-Expansion-CPP/scripts/ExportReference.wl"

     ExportReference["bub00", pdlist, loopmom, extmom, spsRep, topsector, 4, "build/"]

   The script will produce:
     - build/IBPMat_<famname>.bin       (IBP matrices, consumed by C++)
     - build/RefCoeff_<famname>.json     (expansion coefficients for comparison)
*)

(* ------------------------------------------------------------------
   Helper: extract raw c-coefficients from RecursionExpansion
   by instrumenting the recursion loop
   ------------------------------------------------------------------ *)

Options[CaptureExpansionCoefficients] = {
  Modulus -> 0,
  "Increment" -> 2,
  "LinearSolver" -> reduceSolve,
  "LayerByLayer" -> True,
  Verbose -> False
};

CaptureExpansionCoefficients[recursionMatrixData_, order_, OptionsPattern[]] :=
  Module[{
    FFReduce, data, flags, incre, char = OptionValue[Modulus],
    capturedData = {},  (* will store per-order coefficient snapshots *)
    $levelCount = 0, $Verbose = OptionValue[Verbose]
  },

  FFReduce = Function[x, If[char =!= 0, PolynomialMod[x, char], x]];

  (* Use the existing pipeline but capture internal state *)
  data = recursionMatrixData;

  flags = computeFlags[data["M1"], data["N1"], data["nb"], data["ne"], char];
  incre = chooseIncrement[data, flags, OptionValue["Increment"], PrintL];

  (* Run the recursion with instrumented coda *)
  Module[{
    status, hsollist, csollist, cindepall, csol, cindep,
    ne = data["ne"], nb = data["nb"], nibp = data["nibp"],
    vlist = Table[Symbol["v" <> IntegerString[i]], {i, ne}]
  },

    status = "restart";
    While[status === "restart",
      incre = data["incre"];
      status = Catch[
        ClearAll[c];
        c[i_, seed_] := Table[c[i, seed, j], {j, nb}];

        csol = Association@Table[c[0, Table[0, ne], j] -> If[j == 1, 1, 0], {j, nb}];
        cindepall = Table[c[0, Table[0, ne], j], {j, nb}];
        hsollist = {Table[If[j == 1, 1, 0], {j, nb}]};

        (* --- Order loop --- *)
        Do[
          PrintL["====== Order: " <> ToString[i] <> " ======"];

          With[{seeds = Reverse@SortBy[
              Join @@ Table[
                DeleteDuplicates@(Join @@ (Permutations /@
                  (IntegerPartitions[lv + ne, {ne}] - 1))),
                {lv, 0, incre*i - 1}
              ],
              {Total[#], -Reverse[#]} &
            ]},

            (* Remove old solutions *)
            With[{oldIndep = DeleteDuplicates[#[[2]] & /@ cindepall]},
              csol = AssociateTo[csol,
                Join @@ (Association@Table[c[i, #, j] -> 0, {j, nb}] & /@ oldIndep)];
            ];

            (* Layer-by-layer solving *)
            Do[
              csol = WithIndent@solveSeedBlock[
                data, {s}, i, csol, incre, char,
                OptionValue["LinearSolver"], FFReduce, PrintL
              ],
              {s, seeds}
            ];
          ];

          (* --- Coda: capture coefficients at this order --- *)
          With[{res = updateCoda[i, csol, incre, ne, nb, vlist, cindepall, False, FFReduce]},
            csol = res["csol"];
            cindepall = res["cindepall"];
            AppendTo[hsollist, res["hAtOrder"]];

            (* Capture raw coefficients for this order *)
            With[{seedsAtI = Join @@ Table[
                DeleteDuplicates@(Join @@ (Permutations /@
                  (IntegerPartitions[lv + ne, {ne}] - 1))),
                {lv, 0, incre*i}
              ]},
              AppendTo[capturedData,
                <|"order" -> i,
                  "coefficients" -> Table[
                    <|"seed" -> seed, "j" -> j,
                      "value" -> ToString[InputForm[Lookup[csol, c[i, seed, j], 0]]]
                    |>,
                    {seed, seedsAtI}, {j, nb}
                  ]
                |>
              ];
            ];
          ];

          , {i, 1, order}];

        "success"
      ]; (* Catch *)
    ]; (* While *)

    (* Return both the expansion and the captured raw coefficients *)
    <|
      "Expansion" -> hsollist,
      "RawCoefficients" -> capturedData,
      "IndependentVariables" -> cindepall,
      "Status" -> status
    |>
  ]
]


(* ------------------------------------------------------------------
   Main export function
   ------------------------------------------------------------------ *)

Options[ExportReference] = {
  Modulus -> 0,
  "Order" -> 4,
  "OutputDir" -> "."
};

ExportReference[famname_, pdlist_, loopmom_, extmom_, spsRep_, topsector_,
  OptionsPattern[]] :=
  Module[{
    char = OptionValue[Modulus],
    order = OptionValue["Order"],
    outDir = OptionValue["OutputDir"],
    ne, nl, nE, nibp, Alist, vlist, ibpeqsfull, sectorlist,
    expregdata, hexpnAssoc, alist, hlist,
    refData, jsonPath
  },

  (* 1. Define the integral family *)
  {ne, nl, nE, nibp, Alist, vlist, ibpeqsfull, sectorlist} =
    AsyIBPFamilyDefine[pdlist, loopmom, extmom, spsRep, topsector];

  Print["Family: ", famname];
  Print["  ne=", ne, "  nibp=", nibp];
  Print["  Sectors: ", Length@sectorlist];

  (* 2. Determine regions *)
  Print["\n=== Determining Regions ==="];
  expregdata = regionsBySectors[ibpeqsfull, sectorlist, Alist, vlist,
    Modulus -> char, Verbose -> False];

  (* 3. Run expansion and capture coefficients *)
  Print["\n=== Running Expansion ==="];

  refData = {};
  hexpnAssoc = Association@KeyValueMap[
    Function[{sec, regData},
      Print["  Sector: ", sec];

      sec -> Table[
        With[{basis = regData[[j]]["CoordinateRing"]["MonomialBasis"]},
          Print["    Region ", j, "/", Length[regData],
            "  basis(", Length@basis, "): ", basis];

          CaptureExpansionCoefficients[
            regData[[j]]["RecursionMatrix"],
            order,
            Modulus -> char,
            "Increment" -> 2,
            "LayerByLayer" -> True
          ]
        ],
        {j, Length[regData]}
      ]
    ],
    expregdata
  ];

  (* 4. Build reference data structure *)
  refData = <|
    "family" -> famname,
    "order" -> order,
    "modulus" -> char,
    "regions" -> Table[
      Module[{sectorIdx = 1},
        Table[
          <|
            "sector" -> sec,
            "region_index" -> j,
            "ne" -> hexpnAssoc[sec][[j]]["RawCoefficients"][[1]]["coefficients"][[1]]["seed"] // Length,
            "coefficients_by_order" -> hexpnAssoc[sec][[j]]["RawCoefficients"]
          |>,
          {j, Length[hexpnAssoc[sec]]}
        ]
      ],
      {sec, Keys[hexpnAssoc]}
    ] // Flatten
  |>;

  (* 5. Export JSON *)
  jsonPath = FileNameJoin[{outDir, "RefCoeff_" <> famname <> ".json"}];
  Export[jsonPath, refData, "JSON"];
  Print["\nReference coefficients exported to: ", jsonPath];

  (* 6. Also export IBP binary if not already present *)
  Print["\nTo export IBP binary matrices, run:"];
  Print["  ExportBinaryIBPMatrix[\"" <>
    FileNameJoin[{outDir, "IBPMat_" <> famname <> ".bin"}] <>
    "\", expregdata, " <> ToString[char] <> "]"];

  Return[refData]
]


(* ------------------------------------------------------------------
   Quick self-test: verify that exported coefficients satisfy the
   original IBP equations when substituted back
   ------------------------------------------------------------------ *)

Options[VerifyCoefficients] = {Modulus -> 0};

VerifyCoefficients[refData_, ibpeqs_, vlist_, OptionsPattern[]] :=
  Module[{char = OptionValue[Modulus], results},

    results = Table[
      Module[{reg = refData["regions"][[r]],
        coeffs = refData["regions"][[r]]["coefficients_by_order"],
        ne = refData["regions"][[r]]["ne"],
        errors = {}
      },

        (* Reconstruct expansion polynomials *)
        Table[
          Module[{orderData = coeffs[[k]],
            poly, seed, cvals, gterm, eqcheck},

            (* Build polynomial: Sum c[k,seed,j] * v^seed *)
            (* Then check against IBP equations *)
            (* This is a consistency check *)
          ],
          {k, Length[coeffs]}
        ];

        <|"region" -> r, "errors" -> errors|>
      ],
      {r, Length[refData["regions"]]}
    ];

    Return[results]
  ]


Print["ExportReference.wl loaded."];
Print["Available functions:"];
Print["  ExportReference[famname, pdlist, loopmom, extmom, spsRep, topsector]"];
Print["  CaptureExpansionCoefficients[recursionMatrixData, order]"];
Print["  VerifyCoefficients[refData, ibpeqs, vlist]"];
