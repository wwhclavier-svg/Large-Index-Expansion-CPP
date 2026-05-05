(* ::Package:: *)
(* VerifyExpand-MMAExpand.wl — MMA series expansion *)
(* Usage: wolframscript -file VerifyExpand-MMAExpand.wl <famname> [order] *)
(*
    Requires VerifyExpand-Prepare.wl to have been run first (generates .bin + checkpoint).
    If checkpoint is missing, falls back to re-running Prepare steps.
    Output: verify/<fam>/VerifyExpansion-MMAExpansion.m
*)

SetDirectory[DirectoryName[$InputFileName]];
Get["LIEWorkflow.wl"];
Get["VerifyExpand-Prepare.wl"];  (* loads shared utilities *)

(* ---- Parse CLI ---- *)
famname = parseFamilyArg[];
order = If[Length[$ScriptCommandLine] >= 3,
    ToExpression[$ScriptCommandLine[[3]]],
    4
];

config = loadFamilyConfig[famname];
target = targetDir[famname];

Print["========================================"];
Print["VerifyExpand-MMAExpand     fam = ", famname];
Print["  Description: ", config["Description"]];
Print["  Modulus: ", config["Modulus"]];
Print["  Order: ", order];
Print["========================================"];

(* ---- Load checkpoint or re-run Prepare ---- *)
checkpointFile = target <> "PrepareCheckpoint-" <> famname <> ".wdx";
If[FileExistsQ[checkpointFile],
    Print["\n[Step 1] Loading checkpoint: ", checkpointFile];
    data = Get[checkpointFile];
    If[data === Null || !AssociationQ[data] || !KeyExistsQ[data, "Family"],
        Print["  Checkpoint invalid (", Head[data], "), re-running Prepare..."];
        data = Null;
    ,
        Print["  Status: ", data["Status"]];
    ];
];
If[data === Null,
    Print["\n[Step 1b] Running LIEDefineFamily + LIESolveRegions..."];
    data = LIEDefineFamily[
        config["Propagators"], config["LoopMomenta"],
        config["ExternalMomenta"], config["KinematicRules"],
        config["TopSector"],
        "Numeric" -> config["Numeric"],
        Modulus -> config["Modulus"],
        Verbose -> False
    ];
    Print["  NE = ", data["Family", "NE"]];
    Print["  NIBP = ", data["Family", "NIBP"]];

    data = LIESolveRegions[data, Verbose -> False];
    Print["  Regions solved."];
];

(* ---- Expand series ---- *)
Print["\n[Step 2] Expanding series (order=", order, ", LayerByLayer=True)..."];
timeExpand = AbsoluteTiming[
    data = LIEExpandSeries[data,
        "Order" -> order,
        Modulus -> config["Modulus"],
        "Increment" -> 2,
        "LayerByLayer" -> True,
        Verbose -> False
    ];
][[1]];
Print["  Time (Expansion): ", timeExpand, " s"];

hlist = data["Expansion", "HList"];
Print["  HList shape: ", Dimensions[hlist]];

(* ---- Export MMA expansion result ---- *)
Print["\n[Step 3] Exporting MMA expansion..."];

(* Use the top-sector matrix (last one) for SeriesVerify *)
topHList = hlist[[-1, 1]];  (* last matrix, first solution *)
topNB = Length[topHList[[1]]];  (* NB from h^(k=0) component count *)

(* Extract ring metadata from ARegList (1:1 with HList matrices) *)
arList = data["Expansion", "ARegList"];
topCR = arList[[-1]];  (* last matrix = top sector *)

topSectors = Sort[Keys[data[["Regions"]]]];
topSector = Last[topSectors];

qrAFree = topCR["VarIndep"];
qrMinPoly = topCR["MinPoly"];
qrBasis = topCR["MonomialBasis"];
qrVarRule = topCR["VarRule"];

expResult = {{
    <| "Kmax" -> order, "Incre" -> 2,
       "NE" -> data["Family", "NE"], "NB" -> topNB,
       "Nimax" -> data["Family", "NE"] + 4,
       "TopSector" -> topSector,
       "QRingAFree" -> qrAFree,
       "QRingMinPoly" -> qrMinPoly,
       "QRingBasis" -> qrBasis,
       "QRingVarRule" -> qrVarRule,
       "Modulus" -> config["Modulus"],
       "Solutions" -> {
         <| "i" -> 0, "H" -> topHList |>
       }
    |>
}};

mmaFile = target <> "VerifyExpansion-MMAExpansion.m";
Export[mmaFile, expResult, "Text"];
Print["  Top sector: ", topSector, "  NB=", topNB];
Print["  Saved to ", mmaFile];

(* Export timing *)
mmaExpandTiming = <| "Expansion" -> timeExpand |>;
Export[target <> "Compare-MMAExpandTiming-" <> famname <> ".m", mmaExpandTiming, "Text"];

Print["\n========================================"];
Print["VerifyExpand-MMAExpand complete!"];
Print["  Next: wolframscript -file VerifyExpand-Compare.wl ", famname];
Print["========================================"];
