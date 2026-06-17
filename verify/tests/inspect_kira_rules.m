(* inspect_kira_rules.m *)
Get["KiraRuleLoader.wl"];
{kiraRules, kiraReduce} = LoadKiraRules[
    "/root/kira_tests/bub00_new/bub00",
    "d" -> 1/3,
    "Modulus" -> 179424673,
    "NProp" -> 2
];

Print["kiraRules type: ", Head[kiraRules]];
Print["kiraRules length: ", Length[kiraRules]];
Print["kiraRules sample: ", InputForm[kiraRules[[;;5]]]//Short];
Print["kiraRules full: ", InputForm[kiraRules]//Short];
