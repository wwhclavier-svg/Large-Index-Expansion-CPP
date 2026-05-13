(* inspect_fraction_rule.wl *)
$LIECPPPath = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
$VerifyUtilityPath = $LIECPPPath <> "/verify/VerifyUtility/";
$FamilyPath = $LIECPPPath <> "/families/";

SetDirectory[$VerifyUtilityPath];
Get["LIEWorkflow.wl"];
Get["SingularCoordinateRing.wl"];
SetDirectory[$FamilyPath];
Get["FamilyDatabase.wl"];
SetDirectory[$LIECPPPath <> "/verify/DP323_recompute/"];

char = 179424673;
config = $FamilyDatabase["Box"];
workflow = LIEDefineFamily[
  config["Propagators"], config["LoopMomenta"], config["ExternalMomenta"],
  config["KinematicRules"], config["TopSector"],
  "Numeric" -> config["Numeric"], Modulus -> char
];
ibpeqs = workflow["Family", "IBPEqs"];
Alist = workflow["Family", "AList"];
vlist = workflow["Family", "VList"];
sector = {1, 1, 1, 1};

res = regionsBySectors[ibpeqs, {sector}, Alist, vlist,
  Modulus -> char, Verbose -> False];
regions = res[sector];

Do[
  ring = regions[[i]]["CoordinateRing"];
  fr = ring["FractionRule"];
  Print["Region ", i, " FractionRule keys: ", Length[fr], 
    ", Missing count: ", Count[fr, _Missing]];
  If[Count[fr, _Missing] > 0,
    missingKeys = Select[Keys[fr], Head[fr[#]] === Missing &];
    Print["  Missing keys: ", missingKeys];
  ];
, {i, Length[regions]}];
