(* ::Package:: *)
(* Compare-Expand.wl — Generic MMA expansion script *)
(* Usage: wolframscript -file Compare-Expand.wl <famname> [order] *)
(* Output: Compare-MMAResult-<famname>.m *)
(* Reads config from verify/FamilyDatabase/FamilyDatabase.wl *)

SetDirectory[DirectoryName[$InputFileName]];
Get["LIEWorkflow.wl"];
Get[ParentDirectory[DirectoryName[$InputFileName]] <> "/FamilyDatabase/FamilyDatabase.wl"];

(* Parse CLI arguments *)
If[Length[$ScriptCommandLine] < 2,
    Print["ERROR: Missing family name argument."];
    Print["Usage: wolframscript -file Compare-Expand.wl <famname> [order]"];
    PrintAllFamilies[];
    Exit[1]
];
famname = $ScriptCommandLine[[2]];
order = If[Length[$ScriptCommandLine] >= 3,
    ToExpression[$ScriptCommandLine[[3]]],
    4
];

(* Check if family exists *)
If[!KeyExistsQ[$FamilyDatabase, famname],
    Print["ERROR: Unknown family '", famname, "'."];
    PrintAllFamilies[];
    Exit[1]
];

config = $FamilyDatabase[famname];

targetDir = "../" <> famname <> "/";
CreateDirectory[targetDir, CreateIntermediateDirectories -> True];

Print["========================================"];
Print["Compare-Expand     fam = ", famname];
Print["  Description: ", config["Description"]];
Print["  Numeric: ", config["Numeric"]];
Print["  Modulus: ", config["Modulus"]];
Print["  Order: ", order];
Print["========================================"];

(* Step 1: Define Family *)
Print["\n[Step 1] Defining IBP family..."];
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

(* Step 2: Solve Regions *)
Print["\n[Step 2] Solving regions..."];
data = LIESolveRegions[data, Verbose -> False];
Print["  Regions solved."];

(* Step 3: Expand Series *)
Print["\n[Step 3] Expanding series (order=", order, ", LayerByLayer=True)..."];
timeExpand = AbsoluteTiming[
    data = LIEExpandSeries[data,
        "Order" -> order,
        Modulus -> config["Modulus"],
        "Increment" -> 2,
        "LayerByLayer" -> True,
        Verbose -> False
    ];
][[1]];
Print["  Expansion completed."];
Print["  Time (Expansion): ", timeExpand, " s"];

hlist = data["Expansion", "HList"];
Print["  HList shape: ", Dimensions[hlist]];

(* Step 4: Export *)
Print["\n[Step 4] Exporting to Compare-MMAResult-", famname, ".m..."];

expResult = {
  { (* Matrix 1 *)
    <| "Kmax" -> order, "Incre" -> 2,
       "NE" -> data["Family","NE"], "NB" -> 1,
       "Nimax" -> data["Family","NE"] + 2,
       "Solutions" -> {
         <| "i" -> 0, "H" -> hlist[[1, 1]] |>
       }
    |>
  }
};

mmaFile = targetDir <> "Compare-MMAResult-" <> famname <> ".m";
Export[mmaFile, expResult, "Text"];
Print["  Saved to ", mmaFile];

(* Export expansion timing for verification log *)
mmaExpandTiming = <| "Expansion" -> timeExpand |>;
Export[targetDir <> "Compare-MMAExpandTiming-" <> famname <> ".m", mmaExpandTiming, "Text"];

Print["\n========================================"];
Print["Compare-Expand complete!"];
Print["  Next: wolframscript -file Compare-Results.wl ", famname];
Print["========================================"];
