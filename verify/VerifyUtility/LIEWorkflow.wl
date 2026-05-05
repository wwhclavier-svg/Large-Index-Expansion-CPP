(* ::Package:: *)

(* Large Index Expansion - Unified Workflow Package *)
(* Complete workflow: FamilyDefine -> Regions -> Expansion -> Relations *)


BeginPackage["LIEWorkflow`", {
  "LIEUtility`",
  "LIEFamilyDefine`", 
  "LIECoreAlgebra`",
  "LIERegions`",
  "LIEExpand`",
  "LIEReconstruct`"
}]

(* ============================================================ *)
(* Workflow Configuration *)
(* ============================================================ *)

LIEWorkflowConfig::usage = "Association holding workflow configuration.\n"<>
  "Keys: 'Modulus', 'Verbose', 'EnableFieldExtension', 'ExpansionOrder', 'LinearSolver', 'Increment', 'LayerByLayer'"

(* ============================================================ *)
(* Step 1: Family Definition *)
(* ============================================================ *)

LIEDefineFamily::usage = "LIEDefineFamily[pdlist, loopmom, extmom, spsRep, topsector] defines IBP family.\n"<>
  "Returns LIEWorkflowData association with family info."

(* ============================================================ *)
(* Step 2: Region Solving *)
(* ============================================================ *)

LIESolveRegions::usage = "LIESolveRegions[workflowData, opts] or LIESolveRegions[ibpeqs, sectorlist, Alist, vlist, opts] computes expansion regions.\n"<>
  "Returns updated workflowData with region information."

(* ============================================================ *)
(* Step 3: Expansion Reconstruction *)
(* ============================================================ *)

LIEExpandSeries::usage = "LIEExpandSeries[workflowData, order, opts] reconstructs series expansion.\n"<>
  "Returns updated workflowData with expansion coefficients."

(* ============================================================ *)
(* Step 4: Reduction Relations *)
(* ============================================================ *)

LIEGetRelations::usage = "LIEGetRelations[workflowData, rank, opts] reconstructs reduction relations.\n"<>
  "Returns updated workflowData with reduction relations."

(* ============================================================ *)
(* Complete Workflow *)
(* ============================================================ *)

LIECompleteWorkflow::usage = "LIECompleteWorkflow[familyConfig, opts] runs complete workflow in one call.\n"<>
  "familyConfig = <|'Propagators'->..., 'LoopMomenta'->..., 'ExternalMomenta'->..., 'KinematicRules'->..., 'TopSector'->...|>"

LIEWorkflowStep::usage = "LIEWorkflowStep[step, data] executes a single workflow step.\n"<>
  "step can be 'DefineFamily', 'SolveRegions', 'ExpandSeries', 'GetRelations'."

(* ============================================================ *)
(* Data Utilities *)
(* ============================================================ *)

LIENewWorkflowData::usage = "LIENewWorkflowData[] creates empty workflow data association."

LIESaveWorkflowData::usage = "LIESaveWorkflowData[workflowData, filename] saves workflow state."

LIELoadWorkflowData::usage = "LIELoadWorkflowData[filename] loads workflow state."

LIEValidateWorkflowData::usage = "LIEValidateWorkflowData[workflowData, stage] validates data for given stage."

LIEGetRegionSummary::usage = "LIEGetRegionSummary[workflowData] prints summary of computed regions."

LIEGetExpansionSummary::usage = "LIEGetExpansionSummary[workflowData] prints summary of expansion results."


Begin["`Private`"]

(* ============================================================ *)
(* Default Configuration *)
(* ============================================================ *)

$DefaultWorkflowConfig = <|
  "Modulus" -> 0,
  "Verbose" -> True,
  "EnableFieldExtension" -> True,
  "ExpansionOrder" -> 3,
  "LinearSolver" -> reduceSolve,
  "Increment" -> 2,
  "LayerByLayer" -> True,
  "AnsatzMode" -> "Generic",
  "MaxCoefDeg" -> 1 
|>;

(* ============================================================ *)
(* Workflow Data Structure *)
(* ============================================================ *)

(*
LIEWorkflowData structure:
<|
  "Status" -> "New" | "FamilyDefined" | "RegionsSolved" | "ExpansionDone" | "RelationsDone" | "Error",
  "Config" -> <||>,           (* Workflow configuration *)
  "Family" -> <||>,           (* Family definition data *)
  "Regions" -> <||>,           (* Region data (expregdata) *)
  "Expansion" -> <||>,         (* Expansion results *)
  "Relations" -> <||>,         (* Reduction relations *)
  "Timestamps" -> <||>,        (* Timing information *)
  "Errors" -> {}              (* Error messages *)
|>
*)

(* ============================================================ *)
(* Step 1: Family Definition *)
(* ============================================================ *)

Options[LIEDefineFamily] = {
  Modulus -> 0,
  Verbose -> True,
  "Numeric"->{}
};

LIEDefineFamily[pdlist_, loopmom_, extmom_, spsRep_, topsector_, opts:OptionsPattern[]] := Module[
  {ne, nl, nE, nibp, Alist, vlist, ibpeqsfull, sectorlist, timer, config, numeric=OptionValue["Numeric"]},
  
  If[OptionValue[Verbose],
    Print["[LIE-Workflow] Step 1: Defining IBP family..."];
    Print["  - Propagators: ", Length@pdlist];
    Print["  - Loops: ", Length@loopmom];
    Print["  - External legs: ", Length@extmom];
  ];
  
  timer = AbsoluteTiming[
    {ne, nl, nE, nibp, Alist, vlist, ibpeqsfull, sectorlist} = 
      AsyIBPFamilyDefine[pdlist, loopmom, extmom, spsRep, topsector];
  ][[1]];
  
  config = Merge[{$DefaultWorkflowConfig, <|
    "Modulus" -> OptionValue[Modulus],
    "Verbose" -> OptionValue[Verbose]
  |>}, Last];
  
  If[OptionValue[Verbose],
    Print["[LIE-Workflow] Family defined in ", timer, "s"];
    Print["  - Total IBP relations: ", nibp];
    Print["  - Sectors to analyze: ", Length@sectorlist];
  ];
  
  <|
    "Status" -> "FamilyDefined",
    "Config" -> config,
    "Family" -> <|
      "Propagators" -> pdlist,
      "LoopMomenta" -> loopmom,
      "ExternalMomenta" -> extmom,
      "KinematicRules" -> spsRep,
      "TopSector" -> topsector,
      "NE" -> ne,
      "NL" -> nl,
      "NEExt" -> nE,
      "NIBP" -> nibp,
      "AList" -> Alist,
      "VList" -> vlist,
      "IBPEqs" -> ibpeqsfull/.numeric,
      "SectorList" -> sectorlist
    |>,
    "Timestamps" -> <|"FamilyDefined" -> DateString[]|>
  |>
];

(* ============================================================ *)
(* Step 2: Region Solving *)
(* ============================================================ *)

Options[LIESolveRegions] = {
  Modulus -> Automatic,
  "EnableFieldExtension" -> Automatic,
  Verbose -> Automatic
};

LIESolveRegions[workflowData_Association, opts:OptionsPattern[]] := Module[
  {config, family, char, enableFE, verbose, ibpeqs, sectorlist, Alist, vlist, timer, regionData},
  
  (* Validate input *)
  If[!MemberQ[{"FamilyDefined", "RegionsSolved", "ExpansionDone"}, workflowData["Status"]],
    Return[Append[workflowData, {"Status" -> "Error", "Errors" -> {"Family not defined"}}]];
  ];
  
  config = workflowData["Config"];
  family = workflowData["Family"];
  
  (* Get options *)
  char = If[OptionValue[Modulus] === Automatic, config["Modulus"], OptionValue[Modulus]];
  enableFE = If[OptionValue["EnableFieldExtension"] === Automatic, 
    config["EnableFieldExtension"], OptionValue["EnableFieldExtension"]];
  verbose = If[OptionValue[Verbose] === Automatic, config["Verbose"], OptionValue[Verbose]];
  
  (* Extract family data *)
  ibpeqs = family["IBPEqs"];
  sectorlist = family["SectorList"];
  Alist = family["AList"];
  vlist = family["VList"];
  
  If[verbose,
    Print["[LIE-Workflow] Step 2: Solving expansion regions..."];
    Print["  - Sectors: ", Length@sectorlist];
    Print["  - Modulus: ", char];
  ];
  
  timer = AbsoluteTiming[
    regionData = regionsBySectors[ibpeqs, sectorlist, Alist, vlist, 
      Modulus -> char,
      "EnableFieldExtension" -> enableFE,
      Verbose -> verbose
    ];
  ][[1]];
  
  If[verbose,
    Print["[LIE-Workflow] Regions solved in ", timer, "s"];
    Print["  - Non-trivial sectors: ", Length@Keys@regionData];
    Print["  - Total regions: ", Total[Length /@ Values@regionData]];
  ];
  
  Merge[{
    workflowData,
    <|
      "Status" -> "RegionsSolved",
      "Regions" -> regionData,
      "Timestamps" -> Merge[{workflowData["Timestamps"], <|"RegionsSolved" -> DateString[]|>}, Last]
    |>
  }, Last]
];

(* Direct call without workflow data *)
LIESolveRegions[ibpeqs_, sectorlist_, Alist_, vlist_, opts:OptionsPattern[]] := 
  regionsBySectors[ibpeqs, sectorlist, Alist, vlist, opts];

(* ============================================================ *)
(* Step 3: Expansion Reconstruction *)
(* ============================================================ *)

Options[LIEExpandSeries] = {
  "Order" -> Automatic,
  Modulus -> Automatic,
  "LinearSolver" -> Automatic,
  "Increment" -> Automatic,
  "LayerByLayer" -> Automatic,
  Verbose -> Automatic
};

LIEExpandSeries[workflowData_Association, opts:OptionsPattern[]] := Module[
  {config, family, regions, order, char, solver, incre, lbl, verbose, 
   ibpeqs, regionData, timer, hlist, alist},
  
  (* Validate input *)
  If[workflowData["Status"] =!= "RegionsSolved",
    Return[Append[workflowData, {"Status" -> "Error", "Errors" -> {"Regions not solved"}}]];
  ];
  
  config = workflowData["Config"];
  family = workflowData["Family"];
  regions = workflowData["Regions"];
  
  (* Get options *)
  order = If[OptionValue["Order"] === Automatic, config["ExpansionOrder"], OptionValue["Order"]];
  char = If[OptionValue[Modulus] === Automatic, config["Modulus"], OptionValue[Modulus]];
  solver = If[OptionValue["LinearSolver"] === Automatic, config["LinearSolver"], OptionValue["LinearSolver"]];
  incre = If[OptionValue["Increment"] === Automatic, config["Increment"], OptionValue["Increment"]];
  lbl = If[OptionValue["LayerByLayer"] === Automatic, config["LayerByLayer"], OptionValue["LayerByLayer"]];
  verbose = If[OptionValue[Verbose] === Automatic, config["Verbose"], OptionValue[Verbose]];
  
  ibpeqs = family["IBPEqs"];
  regionData = regions;
  
  If[verbose,
    Print["[LIE-Workflow] Step 3: Reconstructing series expansion..."];
    Print["  - Order: ", order];
    Print["  - Regions: ", Total[Length /@ Values@regionData]];
  ];
  
  timer = AbsoluteTiming[
    {hlist, alist} = LIEExpand[
      order, ibpeqs, regionData,
      Modulus -> char,
      "LinearSolver" -> solver,
      "Increment" -> incre,
      "LayerByLayer" -> lbl,
      Verbose -> verbose
    ];
  ][[1]];
  
  If[verbose,
    Print["[LIE-Workflow] Expansion completed in ", timer, "s"];
    Print["  - Coefficients: ", Length@hlist];
  ];
  
  Merge[{
    workflowData,
    <|
      "Status" -> "ExpansionDone",
      "Expansion" -> <|"HList" -> hlist, "ARegList" -> alist, "Order" -> order|>,
      "Timestamps" -> Merge[{workflowData["Timestamps"], <|"ExpansionDone" -> DateString[]|>}, Last]
    |>
  }, Last]
];

(* ============================================================ *)
(* Step 4: Reduction Relations *)
(* ============================================================ *)

Options[LIEGetRelations] = {
  "Rank" -> Automatic,
  Modulus -> Automatic,
  "AnsatzMode" -> Automatic,
  "MaxCoefDeg" -> Automatic,
  Verbose -> Automatic
};

LIEGetRelations[workflowData_Association, opts:OptionsPattern[]] := Module[
  {config, family, expansion, rank, char, ansatzMode, maxDeg, verbose,
   order, hlist, alist, topsec, vlist, timer, relations},
  
  (* Validate input *)
  If[workflowData["Status"] =!= "ExpansionDone",
    Return[Append[workflowData, {"Status" -> "Error", "Errors" -> {"Expansion not done"}}]];
  ];
  
  config = workflowData["Config"];
  family = workflowData["Family"];
  expansion = workflowData["Expansion"];
  
  (* Get options *)
  rank = If[OptionValue["Rank"] === Automatic, config["ExpansionOrder"], OptionValue["Rank"]];
  char = If[OptionValue[Modulus] === Automatic, config["Modulus"], OptionValue[Modulus]];
  ansatzMode = If[OptionValue["AnsatzMode"] === Automatic, config["AnsatzMode"], OptionValue["AnsatzMode"]];
  maxDeg = If[OptionValue["MaxCoefDeg"] === Automatic, config["MaxCoefDeg"], OptionValue["MaxCoefDeg"]];
  verbose = If[OptionValue[Verbose] === Automatic, config["Verbose"], OptionValue[Verbose]];
  
  (* Extract data *)
  order = expansion["Order"];
  hlist = expansion["HList"];
  alist = expansion["ARegList"];
  topsec = family["TopSector"];
  vlist = family["VList"];
  
  If[verbose,
    Print["[LIE-Workflow] Step 4: Reconstructing reduction relations..."];
    Print["  - Rank level: ", rank];
    Print["  - Max coefficient degree: ", maxDeg];
  ];
  
  (* Check if LIEReconstruct is available *)
  If[!MemberQ[Names["*"], "LIEReconstruct"],
    If[verbose,
      Print["[LIE-Workflow] Warning: LIEReconstruct not found. Skipping..."];
      Print["  Please load LIEReconstruct.wl to use this step."];
    ];
    Return[Append[workflowData, {"Status" -> "RelationsDone", "Relations" -> <|"Skipped" -> True|>}]];
  ];
  
  timer = AbsoluteTiming[
    relations = LIEReconstruct[
      rank, maxDeg, order, hlist, alist, topsec, vlist,
      Modulus -> char,
      "AnsatzMode" -> ansatzMode,
      Verbose -> verbose
    ];
  ][[1]];
  
  If[verbose,
    Print["[LIE-Workflow] Relations computed in ", timer, "s"];
    Print["  - Relations found: ", Length@relations];
  ];
  
  Merge[{
    workflowData,
    <|
      "Status" -> "RelationsDone",
      "Relations" -> <|"Relations" -> relations, "Rank" -> rank|>,
      "Timestamps" -> Merge[{workflowData["Timestamps"], <|"RelationsDone" -> DateString[]|>}, Last]
    |>
  }, Last]
];

(* ============================================================ *)
(* Complete Workflow *)
(* ============================================================ *)

Options[LIECompleteWorkflow] = {
  Modulus -> 0,
  Verbose -> True,
  "EnableFieldExtension" -> True,
  "ExpansionOrder" -> 3,
  "ComputeRelations" -> True,
  "LinearSolver" -> reduceSolve,
  "Increment" -> 2,
  "LayerByLayer" -> True
};

LIECompleteWorkflow[familyConfig_Association, opts:OptionsPattern[]] := Module[
  {workflowData, pdlist, loopmom, extmom, spsRep, topsector, computeRel},
  
  (* Extract family configuration *)
  pdlist = familyConfig["Propagators"];
  loopmom = familyConfig["LoopMomenta"];
  extmom = familyConfig["ExternalMomenta"];
  spsRep = familyConfig["KinematicRules"];
  topsector = familyConfig["TopSector"];
  computeRel = OptionValue["ComputeRelations"];
  
  (* Step 1: Define Family *)
  workflowData = LIEDefineFamily[pdlist, loopmom, extmom, spsRep, topsector, 
    "Numeric" -> familyConfig["Numeric"],
    FilterRules[{opts}, Options[LIEDefineFamily]]
  ];
  
  If[workflowData["Status"] === "Error", Return[workflowData]];
  
  (* Step 2: Solve Regions *)
  workflowData = LIESolveRegions[workflowData, 
    FilterRules[{opts}, Options[LIESolveRegions]]
  ];
  
  If[workflowData["Status"] === "Error", Return[workflowData]];
  
  (* Step 3: Expand Series *)
  workflowData = LIEExpandSeries[workflowData, 
    FilterRules[{opts}, Options[LIEExpandSeries]]
  ];
  
  If[workflowData["Status"] === "Error", Return[workflowData]];
  
  (* Step 4: Get Relations (optional) *)
  If[computeRel,
    workflowData = LIEGetRelations[workflowData,
      FilterRules[{opts}, Options[LIEGetRelations]]
    ];
  ];
  
  workflowData
];

(* ============================================================ *)
(* Step-by-Step Interface *)
(* ============================================================ *)

LIEWorkflowStep[step_String, workflowData_Association, opts___] := Switch[step,
  "DefineFamily", LIEDefineFamily[##]& @@ opts,  (* Need special handling *)
  "SolveRegions", LIESolveRegions[workflowData, opts],
  "ExpandSeries", LIEExpandSeries[workflowData, opts],
  "GetRelations", LIEGetRelations[workflowData, opts],
  _, Append[workflowData, {"Status" -> "Error", "Errors" -> {"Unknown step: " <> step}}]
];

(* ============================================================ *)
(* Data Utilities *)
(* ============================================================ *)

LIENewWorkflowData[] := <|
  "Status" -> "New",
  "Config" -> $DefaultWorkflowConfig,
  "Family" -> <||>,
  "Regions" -> <||>,
  "Expansion" -> <||>,
  "Relations" -> <||>,
  "Timestamps" -> <||>,
  "Errors" -> {}
|>;

LIESaveWorkflowData[workflowData_Association, filename_String] := Module[
  {},
  Export[filename, workflowData, "WDX"];
  If[workflowData["Config", "Verbose"],
    Print["[LIE-Workflow] Saved to: ", filename];
  ];
];

LIELoadWorkflowData[filename_String] := Module[
  {data},
  If[!FileExistsQ[filename],
    Print["[LIE-Workflow] Error: File not found: ", filename];
    Return[$Failed];
  ];
  data = Import[filename, "WDX"];
  If[!AssociationQ[data],
    Print["[LIE-Workflow] Error: Invalid data format"];
    Return[$Failed];
  ];
  Print["[LIE-Workflow] Loaded from: ", filename];
  Print["  Status: ", data["Status"]];
  data
];

LIEValidateWorkflowData[workflowData_Association, stage_String] := Module[
  {valid = True, errors = {}},
  
  Switch[stage,
    "DefineFamily",
    If[!KeyExistsQ[workflowData, "Family"],
      valid = False; AppendTo[errors, "Family data missing"];
    ],
    
    "SolveRegions",
    If[workflowData["Status"] =!= "FamilyDefined",
      valid = False; AppendTo[errors, "Family not defined"];
    ],
    
    "ExpandSeries",
    If[workflowData["Status"] =!= "RegionsSolved",
      valid = False; AppendTo[errors, "Regions not solved"];
    ],
    
    "GetRelations",
    If[workflowData["Status"] =!= "ExpansionDone",
      valid = False; AppendTo[errors, "Expansion not done"];
    ],
    
    _,
    valid = False; AppendTo[errors, "Unknown stage: " <> stage];
  ];
  
  <|"Valid" -> valid, "Errors" -> errors|>
];

LIEGetRegionSummary[workflowData_Association] := Module[
  {regions, nSectors, nRegions},
  If[workflowData["Status"] =!= "RegionsSolved" && workflowData["Status"] =!= "ExpansionDone" && workflowData["Status"] =!= "RelationsDone",
    Print["[LIE-Workflow] Regions not computed yet"];
    Return[<||>];
  ];
  
  regions = workflowData["Regions"];
  nSectors = Length@Keys@regions;
  nRegions = Total[Length /@ Values@regions];
  
  Print["=== Region Summary ==="];
  Print["Sectors: ", nSectors];
  Print["Total regions: ", nRegions];
  Print["Regions per sector:"];
  KeyValueMap[
    Print["  ", #1, " -> ", Length@#2, " regions"] &,
    regions
  ];
  
  <|"Sectors" -> nSectors, "TotalRegions" -> nRegions|>
];

LIEGetExpansionSummary[workflowData_Association] := Module[
  {expansion, hlist, order},
  If[workflowData["Status"] =!= "ExpansionDone" && workflowData["Status"] =!= "RelationsDone",
    Print["[LIE-Workflow] Expansion not computed yet"];
    Return[<||>];
  ];
  
  expansion = workflowData["Expansion"];
  hlist = expansion["HList"];
  order = expansion["Order"];
  
  Print["=== Expansion Summary ==="];
  Print["Order: ", order];
  Print["Regions expanded: ", Length@hlist];
  
  <|"Order" -> order, "Regions" -> Length@hlist|>
];


End[] (* Private *)

EndPackage[]
