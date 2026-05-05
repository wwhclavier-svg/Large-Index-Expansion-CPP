(* ::Package:: *)

(* ::Title:: *)
(* M2Kira *)
(* Mathematica Interface for Kira Feynman Integral Reduction *)

(* ::Section:: *)
(* Preamble & Context *)

BeginPackage["M2Kira`"]

$KiraPackageDirectory::usage = "The directory containing the M2Kira package file.";
$KiraExecutablePath::usage = "Path to the Kira executable. Set to a string to override automatic detection.";
$KiraDefaultWorkingDir::usage = "Default working directory for Kira projects. Set to a string path (e.g. \"/root/kira-temp\") or leave as Automatic to use Directory[].";

KiraInputGenerate::usage = "KiraInputGenerate[famName, propagators, loopMom, kinematics, opts] generates all input files (YAMLs and integral lists) required by Kira.";
KiraRun::usage = "KiraRun[famName, opts] executes the Kira binary in the specified working directory.";
KiraRelationImport::usage = "KiraRelationImport[famName, opts] imports Kira reduction results from the generated .m file.";
KiraRelationTest::usage = "KiraRelationTest[rel, sec, vlist, dot, rank, kiraRule, famname, opts] validates Mathematica-generated IBP/LI relations against Kira results.";
(* Internal helper functions kept in Private context to avoid shadowing LargeIndexExpansion.wl *)

(* ::Section:: *)
(* Private *)

Begin["`Private`"]

(* Package directory captured at load time *)
$KiraPackageDirectory = DirectoryName[$InputFileName];
$KiraExecutablePath = Automatic;
$KiraDefaultWorkingDir = Automatic;

(* ------------------------------------------------------------------ *)
(* ::Subsection:: *)
(* Internal Utilities *)
(* ------------------------------------------------------------------ *)

findKira[] := Module[{path, searchPaths},
  (* 1. User-overridden global setting *)
  If[$KiraExecutablePath =!= Automatic && StringQ[$KiraExecutablePath] && FileExistsQ[$KiraExecutablePath],
    Return[$KiraExecutablePath]
  ];
  (* 2. Environment variable KIRAPATH *)
  path = Environment["KIRAPATH"];
  If[path =!= $Failed && FileExistsQ[path], Return[path]];
  (* 3. System PATH via which *)
  path = RunProcess[{"which", "kira"}, "StandardOutput"];
  If[path =!= "" && FileExistsQ[StringTrim[path]], Return[StringTrim[path]]];
  (* 4. Common install / build locations *)
  searchPaths = Flatten[{
    "/usr/local/bin/kira",
    "/usr/bin/kira",
    If[StringQ[$KiraPackageDirectory] && $KiraPackageDirectory =!= "",
      {
        FileNameJoin[{$KiraPackageDirectory, "builddir", "src", "kira", "kira"}],
        FileNameJoin[{$KiraPackageDirectory, "..", "builddir", "src", "kira", "kira"}]
      },
      {}
    ]
  }];
  Do[If[FileExistsQ[p], Return[p]], {p, searchPaths}];
  $Failed
]

yamlValue[expr_] := Module[{str},
  If[NumberQ[expr], Return[ToString[expr, InputForm]]];
  str = ToString[expr, InputForm];
  If[StringMatchQ[str, NumberString], str, "\"" <> str <> "\""]
]

formatPropagator[prop_] := Module[{p = prop},
  If[ListQ[p] && Length[p] == 1, p = p[[1]]];
  If[ListQ[p] && Length[p] == 2,
    {"\"" <> ToString[p[[1]], InputForm] <> "\"", yamlValue[p[[2]]]},
    {"\"" <> ToString[p, InputForm] <> "\"", "0"}
  ]
]

writeFamilyYAML[famName_, propagators_, loopMom_, topSectors_, path_String] := Module[
  {lines, propPairs},
  propPairs = formatPropagator /@ propagators;
  lines = {
    "integralfamilies:",
    "  - name: \"" <> famName <> "\"",
    "    loop_momenta: [" <> StringRiffle[ToString[If[ListQ[#] && Length[#] >= 1, #[[1]], #], InputForm] & /@ loopMom, ", "] <> "]",
    "    top_level_sectors: [" <> StringRiffle[ToString[#] & /@ topSectors, ", "] <> "]",
    "    propagators:"
  };
  lines = Join[lines, ("      - [" <> #[[1]] <> ", " <> #[[2]] <> "]" & /@ propPairs)];
  Export[path, StringRiffle[lines, "\n"] <> "\n", "Text"]
]

writeKinematicsYAML[kin_Association, path_String] := Module[
  {incoming, outgoing, momCons, invariants, scalarRules, sym1, lines},
  incoming = kin["Incoming"] /. _Missing -> {};
  outgoing = kin["Outgoing"] /. _Missing -> {};
  momCons = If[KeyExistsQ[kin, "MomentumConservation"], kin["MomentumConservation"], {}];
  invariants = kin["Invariants"] /. _Missing -> {};
  scalarRules = kin["ScalarRules"] /. _Missing -> {};
  sym1 = kin["SymbolToReplaceByOne"] /. _Missing -> Null;
  
  lines = {
    "kinematics:",
    "  incoming_momenta: [" <> StringRiffle[ToString[#, InputForm] & /@ incoming, ", "] <> "]",
    "  outgoing_momenta: [" <> StringRiffle[ToString[#, InputForm] & /@ outgoing, ", "] <> "]",
    "  momentum_conservation: [" <> StringRiffle[ToString[#, InputForm] & /@ momCons, ", "] <> "]",
    "  kinematic_invariants:"
  };
  lines = Join[lines, ("    - [" <> ToString[#[[1]], InputForm] <> ", " <> ToString[#[[2]]] <> "]" & /@ invariants)];
  lines = Join[lines, {"  scalarproduct_rules:"}];
  lines = Join[lines, ("    - [[" <> ToString[#[[1,1]], InputForm] <> ", " <> ToString[#[[1,2]], InputForm] <> "], " <> yamlValue[#[[2]]] <> "]" & /@ scalarRules)];
  If[sym1 =!= Null,
    lines = Append[lines, "  symbol_to_replace_by_one: " <> ToString[sym1, InputForm]]
  ];
  Export[path, StringRiffle[lines, "\n"] <> "\n", "Text"]
]

writeJobsYAML[famName_, rmax_, smax_, sectorNums_, integralFile_, path_String, extra_, useAutomatic_: False] := Module[
  {lines},
  If[useAutomatic,
    lines = {
      "jobs:",
      " - reduce_sectors:",
      "    select_integrals:",
      "     select_mandatory_list:",
      "      - [" <> famName <> ", " <> integralFile <> "]",
      "    reduce_automatic: true",
      "    run_initiate: true",
      "    run_triangular: true",
      "    run_back_substitution: true",
      " - kira2math:",
      "    target:",
      "     - [" <> famName <> ", " <> integralFile <> "]"
    },
    lines = {
      "jobs:",
      " - reduce_sectors:",
      "    reduce:",
      "     - {topologies: [" <> famName <> "], sectors: [" <> StringRiffle[ToString[#] & /@ sectorNums, ", "] <> "], r: " <> ToString[rmax] <> ", s: " <> ToString[smax] <> "}",
      "    select_integrals:",
      "     select_mandatory_list:",
      "      - [" <> famName <> ", " <> integralFile <> "]",
      "    run_initiate: true",
      "    run_triangular: true",
      "    run_back_substitution: true",
      " - kira2math:",
      "    target:",
      "     - [" <> famName <> ", " <> integralFile <> "]"
    }
  ];
  If[extra =!= {}, lines = Join[lines, extra]];
  Export[path, StringRiffle[lines, "\n"] <> "\n", "Text"]
]

writeIntegralsFile[famName_, seeds_, path_String] := Module[
  {lines},
  lines = (famName <> "[" <> StringRiffle[ToString[#] & /@ #, ", "] <> "]") & /@ seeds;
  Export[path, StringRiffle[lines, "\n"] <> "\n", "Text"]
]

(* ------------------------------------------------------------------ *)
(* ::Subsection:: *)
(* Kira Input File Generation *)
(* ------------------------------------------------------------------ *)

Options[KiraInputGenerate] = {
  "KiraWorkingDir" -> Automatic,
  "RMax" -> 4,
  "SMax" -> 4,
  "Sector" -> Automatic,
  "SelectIntegrals" -> Automatic,
  "IntegralFileName" -> "integrals",
  "TopLevelSectors" -> Automatic,
  "ReduceAutomatic" -> False,
  "ExtraYAML" -> {},
  Verbose -> False
};

KiraInputGenerate[famName_String, propagators_List, loopMom_List, kinematics_Association, opts:OptionsPattern[]] := Module[
  {workDir, famDir, configDir, np, topSectors, sectorVec, sectorNum, seeds, integralPath, jobsPath, familyPath, kinPath},
  
  workDir = OptionValue["KiraWorkingDir"];
  If[workDir === Automatic,
    workDir = If[$KiraDefaultWorkingDir =!= Automatic && StringQ[$KiraDefaultWorkingDir],
      $KiraDefaultWorkingDir,
      Directory[]
    ]
  ];
  If[!DirectoryQ[workDir], CreateDirectory[workDir]];
  famDir = FileNameJoin[{workDir, famName}];
  configDir = FileNameJoin[{famDir, "config"}];
  If[!DirectoryQ[famDir], CreateDirectory[famDir]];
  If[!DirectoryQ[configDir], CreateDirectory[configDir]];
  
  np = Length[propagators];
  topSectors = OptionValue["TopLevelSectors"];
  If[topSectors === Automatic, topSectors = {FromDigits[Reverse[Table[1, {np}]], 2]}];
  
  sectorVec = OptionValue["Sector"];
  If[sectorVec === Automatic, sectorVec = Table[1, {np}]];
  sectorNum = If[ListQ[sectorVec], FromDigits[Reverse[sectorVec], 2], sectorVec];
  If[!ListQ[sectorNum], sectorNum = {sectorNum}];
  
  If[OptionValue["SelectIntegrals"] === Automatic,
    seeds = seedsGenerate[sectorVec, OptionValue["RMax"], OptionValue["SMax"]];
    integralPath = FileNameJoin[{famDir, OptionValue["IntegralFileName"]}];
    writeIntegralsFile[famName, seeds, integralPath];
  ,
    seeds = OptionValue["SelectIntegrals"];
    integralPath = FileNameJoin[{famDir, OptionValue["IntegralFileName"]}];
    writeIntegralsFile[famName, seeds, integralPath];
  ];
  
  familyPath = FileNameJoin[{configDir, "integralfamilies.yaml"}];
  kinPath = FileNameJoin[{configDir, "kinematics.yaml"}];
  jobsPath = FileNameJoin[{famDir, "jobs.yaml"}];
  
  writeFamilyYAML[famName, propagators, loopMom, topSectors, familyPath];
  writeKinematicsYAML[kinematics, kinPath];
  writeJobsYAML[famName, OptionValue["RMax"], OptionValue["SMax"], sectorNum, OptionValue["IntegralFileName"], jobsPath, OptionValue["ExtraYAML"], OptionValue["ReduceAutomatic"]];
  
  If[OptionValue[Verbose] === True, Print["Kira input files generated in: ", famDir]];
  Return[<|
    "FamilyName" -> famName,
    "WorkingDir" -> famDir,
    "IntegralFile" -> OptionValue["IntegralFileName"],
    "JobFile" -> "jobs.yaml",
    "ConfigDir" -> configDir,
    "FamilyPath" -> familyPath,
    "KinematicsPath" -> kinPath,
    "JobsPath" -> jobsPath
  |>];
]

(* ------------------------------------------------------------------ *)
(* ::Subsection:: *)
(* Kira Execution *)
(* ------------------------------------------------------------------ *)

Options[KiraRun] = {
  "KiraExecutable" -> Automatic,
  "FermatPath" -> Automatic,
  "WorkingDir" -> Automatic,
  "JobFile" -> "jobs.yaml",
  "Parallel" -> Automatic,
  Verbose -> True
};

KiraRun[famName_String, opts:OptionsPattern[]] := Module[
  {kiraExe, workDir, env, cmd, result, jobFile, parallel, verbose},
  
  kiraExe = OptionValue["KiraExecutable"];
  If[kiraExe === Automatic, kiraExe = findKira[]];
  If[kiraExe === $Failed, Message[KiraRun::nokira]; Return[$Failed]];
  
  workDir = OptionValue["WorkingDir"];
  If[workDir === Automatic,
    workDir = If[$KiraDefaultWorkingDir =!= Automatic && StringQ[$KiraDefaultWorkingDir],
      FileNameJoin[{$KiraDefaultWorkingDir, famName}],
      FileNameJoin[{Directory[], famName}]
    ]
  ];
  
  jobFile = OptionValue["JobFile"];
  parallel = OptionValue["Parallel"];
  verbose = OptionValue[Verbose];
  
  fermatPath = If[OptionValue["FermatPath"] === Automatic,
    Environment["FERMATPATH"] /. $Failed -> "/usr/local/bin/fermat",
    OptionValue["FermatPath"]
  ];
  
  cmd = {kiraExe, jobFile};
  If[parallel =!= False && parallel =!= None,
    cmd = Append[cmd, If[parallel === Automatic, "--parallel=physical", "--parallel=" <> ToString[parallel]]]
  ];
  
  If[verbose, Print["Running: ", StringRiffle[cmd, " "], " in ", workDir]];
  result = RunProcess[cmd, ProcessDirectory -> workDir, ProcessEnvironment -> <|"FERMATPATH" -> fermatPath|>];
  
  If[result["ExitCode"] =!= 0,
    Message[KiraRun::fail, result["StandardError"]];
    Return[$Failed]
  ];
  
  If[verbose, Print["Kira finished successfully."]];
  Return[result];
]

KiraRun[config_Association?AssociationQ, opts:OptionsPattern[]] := Module[
  {famName, workDir, jobFile, kiraExe, env, cmd, result, parallel, verbose},
  
  famName = Lookup[config, "FamilyName", ""];
  workDir = Lookup[config, "WorkingDir"];
  jobFile = Lookup[config, "JobFile", "jobs.yaml"];
  
  (* opts override config *)
  If[OptionValue["WorkingDir"] =!= Automatic, workDir = OptionValue["WorkingDir"]];
  If[OptionValue["JobFile"] =!= "jobs.yaml", jobFile = OptionValue["JobFile"]];
  
  kiraExe = OptionValue["KiraExecutable"];
  If[kiraExe === Automatic, kiraExe = findKira[]];
  If[kiraExe === $Failed, Message[KiraRun::nokira]; Return[$Failed]];
  
  fermatPath = If[OptionValue["FermatPath"] === Automatic,
    Environment["FERMATPATH"] /. $Failed -> "/usr/local/bin/fermat",
    OptionValue["FermatPath"]
  ];
  
  cmd = {kiraExe, jobFile};
  parallel = OptionValue["Parallel"];
  If[parallel =!= False && parallel =!= None,
    cmd = Append[cmd, If[parallel === Automatic, "--parallel=physical", "--parallel=" <> ToString[parallel]]]
  ];
  
  verbose = OptionValue[Verbose];
  If[verbose, Print["Running: ", StringRiffle[cmd, " "], " in ", workDir]];
  result = RunProcess[cmd, ProcessDirectory -> workDir, ProcessEnvironment -> <|"FERMATPATH" -> fermatPath|>];
  
  If[result["ExitCode"] =!= 0,
    Message[KiraRun::fail, result["StandardError"]];
    Return[$Failed]
  ];
  
  If[verbose, Print["Kira finished successfully."]];
  Return[result];
]

KiraRun::nokira = "Kira executable not found. Please set the path via the 'KiraExecutable' option or add kira to $PATH.";
KiraRun::fail = "Kira exited with an error: `1`.";

(* ------------------------------------------------------------------ *)
(* ::Subsection:: *)
(* Sorting, Seeds, and IBP+LI Generation *)
(* ------------------------------------------------------------------ *)

firstNonZero[list_]:=Module[{i=1},While[list[[i]]===0&&i<=Length@list,i++];
If[i>Length@list,Print["zero list"];Return[-1],Return[i]]];

dotrankseed[top_,d_,r_]:=Module[{pdset,ispset,npd,nisp,dotseeds,rankseeds,seeds},
pdset=Position[#,1][[;;,1]]&@top;ispset=Position[#,0][[;;,1]]&@top;npd=Length@pdset;nisp=Length@ispset;
dotseeds=If[npd==0,{{}},Join@@Table[Join@@(Permutations/@(IntegerPartitions[i+npd,{npd}]-1)),{i,0,d}]];
rankseeds=If[nisp==0,{{}},Join@@Table[Join@@(Permutations/@(IntegerPartitions[i+nisp,{nisp}]-1)),{i,0,r}]];
seeds=Join@@Outer[Normal@SparseArray[Join[Thread[pdset->(1+#1)],Thread[ispset->-#2]],nisp+npd]&,dotseeds,rankseeds,1]
];

MOLabel[index_]:=Join[(List@@index/.{a_/;a>0->1,a_/;a<=0->0})//{Total[#],#}&,{List@@index/.{a_/;a>0:>a-1,a_/;a<=0:>0},List@@index/.{a_/;a>0->0,a_/;a<=0:>-a}}//{Total@#[[1]],Total@#[[2]],#[[2]],#[[1]]}&](* (dot(r),rank(s),s,r) *)

Options[reduceSolve]={Modulus->0};reduceSolve[eqs_,vars_,OptionsPattern[]]:=Module[{char=OptionValue[Modulus],mat,dep,indep,rule,varrev=Reverse@vars},(*Print["eqs: ",Length@eqs,"  vars: ",Length@vars];*)
mat=CoefficientArrays[eqs,varrev]//Join[#[[2]],List/@#[[1]],2]&;
mat=RowReduce[mat,Modulus->char]//DeleteCases[#,a_/;Union[a]=={0}]&(*//Map[rootRationalize,#,{2}]&*);
dep=firstNonZero/@mat;
If[Last@dep>Length@vars,Print["Incompatiable!"];Return[{}]];
indep=Complement[Range[Length@vars],dep];
rule=Thread[varrev[[dep]]->(-mat[[;;,indep]] . varrev[[indep]]-mat[[;;,-1]])];
Return[rule];
];

seedsGenPos[n_Integer,nprop_]:=If[nprop==0,{{}},Join@@Table[DeleteDuplicates@(Join@@(Permutations/@(IntegerPartitions[i+nprop,{nprop}]-1))),{i,0,n}]];
seedsGenPos[n_List,nprop_]:=If[nprop==0,{{}},Join@@Table[DeleteDuplicates@(Join@@(Permutations/@(IntegerPartitions[n[[i]]+nprop,{nprop}]-1))),{i,Length@n}]];

seedsGenerate[sector_List,r_?IntegerQ,s_?IntegerQ]:=Module[{corner,inddots,indrank,seeds,npd,sec},
npd=Length@sector;sec=Flatten@Position[sector,1];
corner=Plus@@(UnitVector[npd,#]&/@sec);
inddots=(# . (UnitVector[npd,#]&/@sec))&/@seedsGenPos[r,Length[sec]];
indrank=(# . (UnitVector[npd,#]&/@Complement[Range[npd],sec]))&/@seedsGenPos[s,npd-Length[sec]];
seeds=Outer[(corner+#1-#2)&,inddots,indrank,1]//Flatten[#,1]&;
Return[seeds];
];

(* ordering of FIs according to Kira *)
orderingfi[fi_]:={Total@#[[1]],#[[1]],Total@#[[2]],Total@#[[3]],-#[[3]],-#[[2]]}&@{#/.{b_?Positive->1,b_?Negative->0},Select[List@@#,Positive],-Select[List@@#,NonPositive]}&@(List@@fi)
BLSort[FIlist_]:=SortBy[FIlist,orderingfi];

mon2F[poly_,zlist_,fsymb_]/;Head[poly]=!=List:=Total[(fsymb@@Exponent[#,zlist])*(#/.(#->1&/@zlist))&/@If[Head[poly]===Plus,List@@poly,List@poly]];
mon2F[vec_List,zlist_,fsymb_]:=mon2F[#,zlist,fsymb]&/@vec;
f2Mon[expr_,zlist_,fsymb_]:=expr/.fsymb[ind__]:>Times@@Thread@Power[zlist,{ind}];

SP2PD[pdlist_,loopmom_,extmom_,spsRep_]:=Module[{scalarproducts,pd2sp,sp2pd,sp2pdrule,npd=Length@pdlist},
scalarproducts=Join[Union@Flatten@Outer[Times,loopmom,loopmom],Flatten@Outer[Times,loopmom,extmom]];
pd2sp=Expand[pdlist]/.Thread[Rule[scalarproducts,Array["sp",npd]]]/.spsRep;(*express pdlist in terms of sp[i]*)
sp2pd=LinearSolve[#[[2]],-#[[1]]+Array["z",npd]]&@CoefficientArrays[Expand[pdlist]/.Thread[Rule[scalarproducts,Array["sp",npd]]]/.spsRep,Array["sp",npd]];(*express scalarproduces in terms of z[i]*)
sp2pdrule=Thread[Rule[scalarproducts,sp2pd]];
Return[sp2pdrule];
];

derivL[pdlist_,loopmom_,extmom_,spsRep_,sprule_]:=(Expand[Expand[Outer[Times,D[#,{loopmom}],Join[loopmom,extmom]]]/.spsRep/.sprule])&/@pdlist
derivP[pdlist_,loopmom_,extmom_,spsRep_,sprule_]:=(Expand[Expand[Outer[Times,If[extmom=!={},D[#,{extmom}],{}],extmom]]/.spsRep/.sprule])&/@pdlist;
(*derivative of propagator denominators in terms of loop momenta or external momenta: q_j dD_i/dl_k, p_j dD_i/dp_k*)

(* Set default option for genIBP so that it works out of the box *)
Options[genIBP] = {"dimension" -> d};

genIBP[derivl_]:=Module[{genRel,dim=OptionValue["dimension"]},
genRel[j_Integer,k_Integer,ind_,derivmat_]:=Module[{eqmon,eqf,npd=Length@derivmat},
eqmon=If[j==k,$d,0]+Sum[-ind[[i]]*derivmat[[i,j,k]]/"z"[i],{i,npd}];
eqmon=Expand[eqmon/(Times@@Thread@Power["z"/@Range[npd],ind])];
eqf=mon2F[eqmon,"z"/@Range[npd],"f"]/.{"f"[vec__]:>"f"@@(-{vec})};
Return[eqf//Collect[#,_f]&];
];
genRel[#[[1]],#[[2]],Array["n"<>ToString[#]&,Length@derivl],derivl]&/@(Join@@Array[{#1,#2}&,Dimensions[derivl][[2;;3]]])
];

genLI[derivp_,extmom_,spsRep_]:=Module[{genRel,ne=Length@extmom},
genRel[j_Integer,k_Integer,ind_,derivmat_]:=Module[{eqLI,eqLIf,npd=Length@derivmat},eqLI=Sum[-ind[[i1]]*Sum[(Times@@extmom[[{j,i2}]]/.spsRep)*derivmat[[i1,i2,k]]-(Times@@extmom[[{k,i2}]]/.spsRep)*derivmat[[i1,i2,j]],{i2,ne}]/"z"[i1],{i1,npd}];
eqLI=Expand[eqLI/(Times@@Thread@Power["z"/@Range[npd],ind])];
eqLIf=mon2F[eqLI,"z"/@Range[npd],"f"]/.{"f"[vec__]:>"f"@@(-{vec})};
Return[eqLIf];
];
genRel[#[[1]],#[[2]],Array["n"<>ToString[#]&,Length@derivp],derivp]&/@(Join@@Table[Table[{i,j},{j,i+1,ne}],{i,ne}])
];

relationGeneral[famname_,pdlist_,loopmom_,extmom_,spsRep_,OptionsPattern[]]:=Module[{sprule,derivl,derivp,relationgeneral,relationspecific,npd=Length@pdlist},
sprule=SP2PD[pdlist,loopmom,extmom,spsRep];
derivl=derivL[pdlist,loopmom,extmom,spsRep,sprule];
derivp=derivP[pdlist,loopmom,extmom,spsRep,sprule];
relationgeneral=genIBP[derivl]~Join~genLI[derivp,extmom,spsRep];
Return[relationgeneral/.{"f"[ind__]:>j@@Join[{famname},{ind}]}];
]
Options[relationSpecific]={Cut->{}};
relationSpecific[famname_,seedslist_,pdlist_,loopmom_,extmom_,spsRep_,OptionsPattern[]]:=Module[{sprule,derivl,derivp,relationgeneral,relationspecific,npd=Length@pdlist,cut=OptionValue[Cut]},
sprule=SP2PD[pdlist,loopmom,extmom,spsRep];
derivl=derivL[pdlist,loopmom,extmom,spsRep,sprule];
derivp=derivP[pdlist,loopmom,extmom,spsRep,sprule];
relationgeneral=genIBP[derivl]~Join~genLI[derivp,extmom,spsRep];
relationspecific=Join@@Table[(relationgeneral/.Thread[Array["n"<>ToString[#]&,npd]->seedslist[[i]]])/.{"f"[ind__]/;Or@@(#<=0&/@{ind}[[cut]])->0},{i,Length@seedslist}];
relationspecific=DeleteCases[relationspecific,0];
Return[relationspecific/.{"f"[ind__]:>j@@Join[{famname},{ind}]}];
]

(* ------------------------------------------------------------------ *)
(* ::Subsection:: *)
(* Result Import *)
(* ------------------------------------------------------------------ *)

Options[kiraRelationImport]={"WorkingDir"->Automatic,"IntegralFile"->"integrals","ExplicitFamilyName"->False};
kiraRelationImport[famname_,OptionsPattern[]]:=Module[{kiraRule,kiraMI,workDir,integralFile,resultPath},
  workDir=OptionValue["WorkingDir"];
  If[workDir===Automatic,
    workDir=If[$KiraDefaultWorkingDir=!=Automatic&&StringQ[$KiraDefaultWorkingDir],
      FileNameJoin[{$KiraDefaultWorkingDir,famname}],
      FileNameJoin[{Directory[],famname}]
    ]
  ];
  integralFile=OptionValue["IntegralFile"];
  resultPath=FileNameJoin[{workDir,"results",famname,"kira_"<>integralFile<>".m"}];
  If[!FileExistsQ[resultPath],
    Message[kiraRelationImport::nofile,resultPath];
    Return[$Failed]
  ];
  kiraRule=Get[resultPath];
  kiraMI=Union@Cases[kiraRule[[;;,2]],a_/;ToString@Head[a]===famname,Infinity]/. ToExpression[famname][a__]:>Global`j[famname,{a}];
  kiraRule=kiraRule/. ToExpression[famname][a__]:>Global`j[famname,{a}];
  Print["MI in kira: ",kiraMI];
  If[OptionValue["ExplicitFamilyName"]==False,{kiraRule,kiraMI}={kiraRule,kiraMI}/.{Global`j[a_,b_]:>Global`j@@b}];
  Return[{kiraRule,kiraMI}];
]
kiraRelationImport::nofile="Kira result file not found: `1`";

(* Public wrapper with capitalised name, delegating to the original implementation *)
Options[KiraRelationImport] = Options[kiraRelationImport];
KiraRelationImport[famname_, opts:OptionsPattern[]] := kiraRelationImport[famname, opts]

KiraRelationImport[config_Association?AssociationQ, opts:OptionsPattern[]] := Module[
  {famname, workDir, integralFile, resultPath, kiraRule, kiraMI},
  famname = Lookup[config, "FamilyName"];
  workDir = Lookup[config, "WorkingDir"];
  integralFile = Lookup[config, "IntegralFile", "integrals"];
  
  (* opts override config *)
  If[OptionValue["WorkingDir"] =!= Automatic, workDir = OptionValue["WorkingDir"]];
  If[OptionValue["IntegralFile"] =!= "integrals", integralFile = OptionValue["IntegralFile"]];
  
  resultPath = FileNameJoin[{workDir, "results", famname, "kira_" <> integralFile <> ".m"}];
  If[!FileExistsQ[resultPath],
    Message[KiraRelationImport::nofile, resultPath];
    Return[$Failed]
  ];
  kiraRule = Get[resultPath];
  kiraMI = Union@Cases[kiraRule[[;;,2]], a_/;ToString@Head[a]===famname, Infinity] /. ToExpression[famname][a__] :> Global`j[famname, {a}];
  kiraRule = kiraRule /. ToExpression[famname][a__] :> Global`j[famname, {a}];
  Print["MI in kira: ", kiraMI];
  If[OptionValue["ExplicitFamilyName"] == False,
    {kiraRule, kiraMI} = {kiraRule, kiraMI} /. {Global`j[a_, b_] :> Global`j @@ b}
  ];
  Return[{kiraRule, kiraMI}];
]

(* ------------------------------------------------------------------ *)
(* ::Subsection:: *)
(* Result Verification *)
(* ------------------------------------------------------------------ *)

Options[kiraRelationTest]={Mode->0,Modulus->0,"RelationDepth"->2,Verbose->False};
kiraRelationTest[rel_,sec_,vlist_,dot_,rank_,kiraRule_,famname_,OptionsPattern[]]:=Module[{depth,testseed,reln,reln1,lmrel,mode=OptionValue[Mode],char=OptionValue[Modulus],printF},depth=Max[(Max[#]-Min[#])&@Cases[#,g[a__]:>Total[{a}-vlist],Infinity]&/@rel];
testseed=dotrankseed[sec,dot,Max[1,rank-depth+1]];
Table[
Print["seed: ",seed];
reln=(rel/.{n->0}/.Thread[vlist->seed]/.{g[a__]/;And@@(#<=0&/@{a}):>0});reln=reln//Select[#,Cases[#,g[a__]/;Total[{a}/.{b_/;b<0:>-b,b_/;b>0:>0}]>rank,Infinity]=={}&]&//DeleteCases[#,0]&;reln=If[reln=!={},
reduceSolve[#==0&/@reln,SortBy[Union@Cases[reln,_g,Infinity],MOLabel[List@@#]&]]//#[[1]]-#[[2]]&/@#&,{}];
(* Drop: intg. in all zero sector (0...0); eqs. with intg higher rank *)
lmrel=If[#=!={},Last[#],#]&@SortBy[Union@Cases[#,_g,{0,Infinity}],MOLabel[List@@#]&]&/@reln;
reln1=reln//Collect[Expand[#/.(kiraRule/.ToExpression[famname][a__]:>g[a])/.{g[a__]/;And@@(#<=0&/@{a})->0}],_g]&;
If[char=!=0,reln1=PolynomialMod[reln1,char]];
(* output *)
Switch[mode,0,
If[reln=!={},
Print[Join[List/@lmrel,List/@reln,2]//MatrixForm,"  ",Factor@reln1//MatrixForm,"\n--------"]],
1,
If[reln=!={}&&Union[reln1]=!={0},
Print[Join[List/@lmrel,List/@reln,2]//MatrixForm,"  ",Factor@reln1//MatrixForm,"\n--------"]]
],{seed,testseed}];
]

(* Public wrapper with capitalised name, delegating to the original implementation *)
KiraRelationTest[rel_, sec_, vlist_, dot_, rank_, kiraRule_, famname_, opts:OptionsPattern[]] := kiraRelationTest[rel, sec, vlist, dot, rank, kiraRule, famname, opts]

(* ------------------------------------------------------------------ *)
(* ::Section:: *)
(* End of Package *)
(* ------------------------------------------------------------------ *)

End[]
EndPackage[]
