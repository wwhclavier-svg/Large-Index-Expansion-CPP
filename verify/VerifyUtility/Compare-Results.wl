(* ::Package:: *)
(* Compare-Results.wl — Generic C++ vs MMA expansion comparison *)
(* Usage: wolframscript -file Compare-Results.wl <famname> *)
(* Auto-searches MMA-Mini, C++ project root, and build/ for result files *)

SetDirectory[DirectoryName[$InputFileName]];

(* Parse CLI argument *)
If[Length[$ScriptCommandLine] < 2,
    Print["ERROR: Missing family name argument."];
    Print["Usage: wolframscript -file Compare-Results.wl <famname>"];
    Exit[1]
];
famname = $ScriptCommandLine[[2]];

modulus = Prime[10000000];

Print["========================================"];
Print["Compare-Results     fam = ", famname];
Print["  Modulus = ", modulus];
Print["========================================"];

(* Search for a file across multiple base directories *)
findFile[baseName_] := Module[{searchPaths, found},
    searchPaths = {"",
                   "../../",
                   "../../build/",
                   "../" <> famname <> "/"};
    found = SelectFirst[searchPaths, FileExistsQ[# <> baseName] &];
    If[MissingQ[found], baseName, found <> baseName]
];

baseCPP = "Compare-CPPResult-" <> famname <> ".m";
baseMMA = "Compare-MMAResult-" <> famname <> ".m";
cppFile = findFile[baseCPP];
mmaFile = findFile[baseMMA];

(* Step 1: Read C++ result *)
Print["\n[Step 1] Reading C++ result from: ", cppFile];
If[!FileExistsQ[cppFile],
    Print["[FAIL] C++ result file not found: ", cppFile];
    Print["  Searched: . , ../../ , ../../build/"];
    Print["  Run: cd ../.. && ./build/test_expandFF ", famname];
    Exit[1];
];
Get[cppFile];
hlistCPP = $ExpansionResults[[1, 1, "Solutions", 1, "H"]];
Print["  C++ result loaded. Orders: ", Length[hlistCPP] - 1];

(* Step 2: Read MMA result *)
Print["\n[Step 2] Reading MMA result from: ", mmaFile];
If[!FileExistsQ[mmaFile],
    Print["[FAIL] MMA result file not found: ", mmaFile];
    Print["  Searched: . , ../../ , ../../build/"];
    Print["  Run: wolframscript -file Compare-Expand.wl ", famname];
    Exit[1];
];
Get[mmaFile];
hlistMMA = $ExpansionResults[[1, 1, "Solutions", 1, "H"]];
Print["  MMA result loaded. Orders: ", Length[hlistMMA] - 1];

(* Step 3: Compare each order *)
Print["\n[Step 3] Comparing coefficients..."];

maxk = Min[Min[Length[hlistCPP], Length[hlistMMA]] - 1, 4];
allMatch = True;
For[k = 0, k <= maxk, k++,
    cppPoly = hlistCPP[[k + 1]];
    mmaPoly = hlistMMA[[k + 1]];

    diff = PolynomialMod[Expand[mmaPoly - cppPoly], modulus];

    If[diff === 0,
        Print["  [MATCH] k=", k],
        Print["  [MISMATCH] k=", k];
        Print["    C++ = ", cppPoly];
        Print["    MMA = ", mmaPoly];
        Print["    diff = ", diff];
        allMatch = False;
    ];
];

If[maxk < 4,
    Print["  (only orders 0..", maxk, " were available in both results)"]
];

Print["\n========================================"];
If[allMatch,
    Print["[PASS] All orders match! C++ and MMA are consistent."],
    Print["[FAIL] Some orders do not match. See details above."]
];
Print["========================================"];
If[!allMatch, Exit[1]];
