(* ::Package:: *)
(* KiraRuleLoader.wl *)
(* 统一加载 Kira 规则，包括显式约化规则、平凡 sector 规则、零指标规则 *)
(* 用法: Get["KiraRuleLoader.wl"]; {kiraRules, kiraReduce} = LoadKiraRules[kiraBaseDir]; *)

BeginPackage["KiraRuleLoader`"]

j::usage = "Integral symbol used in Kira reduction rules.";

LoadKiraRules::usage = "LoadKiraRules[kiraBaseDir, opts] loads all Kira rules (explicit reductions + trivial sector rules + zero-index rule) from a Kira project directory. Returns {allRules, kiraReduceFunc}.";

Begin["`Private`"]

LoadKiraRules[kiraBaseDir_String, OptionsPattern[]] := Module[
  {modulus, dVal, kiraFile, kiraRaw, kiraExplicit, zeroIndexRule,
   trivialFile, trivialSectors, nprop, trivialRules, allRulesSymbolic, allRulesFF, reduceFunc},

  modulus = OptionValue["Modulus"] /. Automatic -> Prime[10000000];
  dVal = OptionValue["d"] /. Automatic -> 1/3;
  nprop = OptionValue["NProp"] /. Automatic -> 2;
  familyName = OptionValue["FamilyName"] /. Automatic -> "bub00";

  (* 1. 加载 Kira 显式约化规则 *)
  kiraFile = FileNameJoin[{kiraBaseDir, "results", familyName, "kira_integrals.m"}];
  If[!FileExistsQ[kiraFile],
    Print["ERROR: Kira rules file not found: ", kiraFile];
    Return[{$Failed, $Failed}]
  ];
  (* 加载 Kira 规则文件，确保 bub00 被正确替换为 j *)
  (* 使用字符串替换避免上下文问题 *)
  kiraText = Import[kiraFile, "Text"];
  kiraText = StringReplace[kiraText, familyName -> "j"];
  kiraExplicit = ToExpression[kiraText];

  (* 2. 精确零指标规则 *)
  (* 用户指定: j[a__] 当所有 a <= 0 时 -> 0 *)
  zeroIndexRule = j[a__]/;!Or@@({a}/.{b_/;b>0->True,b_/;b<=0->False})->0;

  (* 3. Kira 平凡 sector 规则 *)
  (* 从 Kira 的 sectormappings 目录读取 trivial sector 信息 *)
  trivialFile = FileNameJoin[{kiraBaseDir, "sectormappings", familyName, "trivialsector"}];
  trivialRules = {};
  If[FileExistsQ[trivialFile],
    trivialSectors = ToExpression[StringSplit[Import[trivialFile, "Text"], ","]];
    (* 排除 sector 0 (全零，已由零规则覆盖) *)
    trivialSectors = Select[trivialSectors, # > 0 &];
    trivialRules = Flatten[Table[
      Module[{secVec, conds, rule},
        secVec = Reverse[IntegerDigits[secNum, 2, nprop]];
        (* 生成条件: ai > 0 当 si=1, ai <= 0 当 si=0 *)
        (* 对于 nprop=2, sector {1,0} -> a>0 && b<=0 *)
        (* 对于 nprop=2, sector {0,1} -> a<=0 && b>0 *)
        conds = Table[
          If[secVec[[i]] == 1,
            Symbol["a" <> ToString[i]] > 0,
            Symbol["a" <> ToString[i]] <= 0
          ],
          {i, nprop}
        ];
        (* 构建 Mathematica 模式规则 *)
        If[nprop == 2,
          If[secNum == 1,
            j[a_, b_] /; a > 0 && b <= 0 -> 0,
            j[a_, b_] /; a <= 0 && b > 0 -> 0
          ],
          (* 通用 nprop 规则 *)
          Module[{vars = Table[ToExpression["a" <> ToString[i]], {i, nprop}]},
            j[Sequence @@ (Pattern[#, Blank[]] & /@ vars)] /; And @@ Table[
              If[secVec[[i]] == 1, ToExpression["a" <> ToString[i]] > 0, ToExpression["a" <> ToString[i]] <= 0],
              {i, nprop}
            ] -> 0
          ]
        ]
      ],
      {secNum, trivialSectors}
    ]];,
    Print["WARNING: Trivial sector file not found: ", trivialFile];
  ];

  (* 合并所有规则（保持符号形式，含 d） *)
  allRulesSymbolic = Join[kiraExplicit, {zeroIndexRule}, trivialRules];

  (* 返回约化函数 *)
  (* 注意: 不在规则层面取模，而是在应用规则后代入 d 再取模 *)
  (* 这样可以避免 PolynomialMod 破坏有理数系数规则 *)
  (* 使用 With 将 dVal/modulus/allRulesSymbolic 注入 Function 体， *)
  (* 并通过 Set (=) 让 Module 返回纯函数对象而非局部符号，避免 Module 清除 *)
  With[{dv = dVal, mod = modulus, rules = allRulesSymbolic},
    reduceFunc = Function[{expr},
      Module[{result},
        result = FixedPoint[(# /. rules) &, expr, 100];
        (* Kira规则中d位于Global上下文，需显式指定以匹配 *)
        result = result /. Global`d -> dv;
        (* 对有理数系数进行模运算 *)
        result = result /. Rational[p_, q_] :> Mod[p * PowerMod[q, -1, mod], mod];
        result = PolynomialMod[result, mod]
      ]
    ]
  ];

  Print["Kira rules loaded: "];
  Print["  Explicit rules:   ", Length[kiraExplicit]];
  Print["  Zero-index rule:  1"];
  Print["  Trivial sectors:  ", Length[trivialRules]];
  Print["  Total rules:      ", Length[allRulesSymbolic]];

  (* 返回符号规则列表和约化函数 *)
  {allRulesSymbolic, reduceFunc}
];

Options[LoadKiraRules] = {
  "Modulus" -> Automatic,
  "d" -> Automatic,
  "NProp" -> Automatic,
  "FamilyName" -> Automatic
};

End[]
EndPackage[]
