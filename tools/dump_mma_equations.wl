(* Dump MMA IBP equations and A/B equations for comparison with C++ *)
SetDirectory[DirectoryName[$InputFileName] <> "/../verify/VerifyUtility"];

Get["LIEWorkflow.wl"];
Get[ParentDirectory[DirectoryName[$InputFileName] <> "/../verify"] <> "/FamilyDatabase/FamilyDatabase.wl"];

config = $FamilyDatabase["bub00"];

Print["=== MMA bub00 IBP Equations ==="];
Print["Modulus: ", config["Modulus"]];
Print["Numeric: ", config["Numeric"]];

(* Step 1: Define Family *)
data = LIEDefineFamily[
    config["Propagators"], config["LoopMomenta"],
    config["ExternalMomenta"], config["KinematicRules"],
    config["TopSector"],
    "Numeric" -> config["Numeric"],
    Modulus -> config["Modulus"],
    Verbose -> False
];

ne = data["Family", "NE"];
nibp = data["Family", "NIBP"];
ibpeqs = data["Family", "IBPEqs"];
vlist = data["Family", "VList"];
Alist = data["Family", "AList"];

Print["NE=", ne, " NIBP=", nibp];
Print["VList: ", vlist];

Print["--- Raw IBP Equations ---"];
Do[
  Print["Eq[", i, "]: ", InputForm[ibpeqs[[i]]]];
, {i, nibp}];

(* Build A/B equations (same as expRegSolve2 for top sector) *)
AlistCF = Alist /. A[i_] :> "A"[i];
aeqs0 = Coefficient[ibpeqs, "n"] /. "g"[a__] :> (Times @@ Thread @ Power[AlistCF, {a} - vlist]);
aeqs = Join[aeqs0 /. {1/"A"[a_] :> "B"[a]}, Table["A"[i] "B"[i] - 1, {i, ne}]];

Print["--- A/B Equations (top sector, v=all-1) ---"];
Do[
  Print["AB[", i, "]: ", InputForm[aeqs[[i]]]];
, {i, Length[aeqs]}];

(* Compute GB *)
Print["--- Groebner Basis ---"];
agb = GroebnerBasis[aeqs, Join @@ {{"B"[1], "B"[2]}, {"A"[1], "A"[2]}}, Modulus -> config["Modulus"]];
Do[
  Print["GB[", i, "]: ", InputForm[agb[[i]]]];
, {i, Length[agb]}];

(* Compute VarRule *)
Print["--- VarRule ---"];
xdrule = Solve[# == 0 & /@ agb, {"A"[1], "A"[2], "B"[1], "B"[2]}, Modulus -> config["Modulus"]];
Print["Solutions: ", xdrule];
