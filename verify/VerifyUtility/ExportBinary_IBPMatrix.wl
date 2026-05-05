(* ::Package:: *)

BeginPackage["ExportIBPMatrixBinary`"]

ExportBinaryIBPMatrix::usage = "ExportBinaryIBPMatrix[filename, expregdata, mod] exports IBP matrices to binary format."
ExportBinaryRingData::usage = "ExportBinaryRingData[filename, expregdata, Alist, ne, mod] exports ring data (matrices A and Ainv) to binary format."

Begin["`Private`"]

monomialRulesPower[expr_,var_]:=If[var=!={},
Association[Join@@
(If[Head[#]===SparseArray,
Table[Length@Cases[#[[1]],i,Infinity],{i,Length@var}]->#[[2]]&/@Most@ArrayRules[#],{Table[0,Length@var]->#}]&/@CoefficientArrays[expr,var])],<|{}->expr|>];
monomialRulesPower[expr_List,var_]:=Module[{rulelist,monlist},
rulelist=monomialRulesPower[#,var]&/@expr;monlist=Union@@(Keys/@rulelist);
Association@Table[mon->Table[If[MemberQ[Keys@rule,mon],rule[mon],0],{rule,rulelist}],{mon,monlist}]];

(* \:8f85\:52a9\:51fd\:6570\:ff1a\:5c06\:7a00\:758f\:77e9\:9635\:8f6c\:6362\:4e3a\:4e8c\:8fdb\:5236\:6570\:636e\:5217\:8868\:ff08\:4e0d\:5305\:62ec\:5b58\:5728\:6807\:5fd7\:ff09 *)
sparseToBinary[mat_SparseArray] := Module[
    {dims = mat["Dimensions"],
     rowPtr = mat["RowPointers"],
     colIdx = Flatten[mat["ColumnIndices"]],
     vals = mat["NonzeroValues"]},
    {
        Length[dims], dims,
        Length[rowPtr], rowPtr,
        Length[colIdx], colIdx,
        Length[vals], vals
    }
];

(* \:4e3b\:5bfc\:51fa\:51fd\:6570 *)
ExportBinaryIBPMatrix[filename_String, expregdata_, mod_Integer] := Module[
    {recursionMatrixList, stream, numRegs, operatorOrder, regData,
     nibp, ne, nb, incre, exists, opData, binaryData, dimsLen, dims, rowPtrLen, rowPtr, colIdxLen, colIdx, valsLen, vals, typeList},

    (* \:63d0\:53d6\:6240\:6709 recursionMatrix\:ff0c\:5047\:8bbe expregdata \:7ed3\:6784\:4e0e IBPMatrix_Export \:517c\:5bb9 *)
    recursionMatrixList = Cases[Flatten[Values[expregdata]], _Association];
    (* \:6216\:8005\:66f4\:7cbe\:786e\:ff1a\:6bcf\:4e2a\:5143\:7d20\:662f Association\:ff0c\:4e14\:542b\:6709 "RecursionMatrix" \:952e *)
    recursionMatrixList = (#["RecursionMatrix"] &) /@ recursionMatrixList;

    numRegs = Length[recursionMatrixList];
    If[numRegs == 0, Message[ExportBinaryIBPMatrix::nodata]; Return[$Failed]];
    Print["#reg: ",numRegs];

    stream = OpenWrite[filename, BinaryFormat -> True];
    If[stream === $Failed, Return[$Failed]];

    (* \:5199\:5165\:9b54\:6570 "IBP1" *)
    BinaryWrite[stream, "IBP1", "Character8", ByteOrdering -> +1];

    (* \:5199\:5165\:7ec4\:6570 *)
    BinaryWrite[stream, numRegs, "Integer32", ByteOrdering -> +1];

    (* \:7b97\:5b50\:987a\:5e8f *)
    operatorOrder = {"M1", "N1", "K1", "F0", "F2", "K1s", "K2s", "F2s"};

    Do[
        regData = recursionMatrixList[[i]];

        (* \:4ece regData \:63d0\:53d6\:5143\:6570\:636e *)
        (* \:4f18\:5148\:4ece M1 \:83b7\:53d6\:ff0c\:82e5\:6ca1\:6709\:5219\:4ece N1 \:7b49\:83b7\:53d6 *)
		Which[
		    KeyExistsQ[regData, "M1"],
		    nibp = regData["M1"]["Dimensions"][[1]];
		    ne   = regData["M1"]["Dimensions"][[2]];
		    nb   = regData["M1"]["Dimensions"][[3]],
		    KeyExistsQ[regData, "N1"],
		    nibp = regData["N1"]["Dimensions"][[1]];
		    ne   = regData["N1"]["Dimensions"][[2]];
		    nb   = regData["N1"]["Dimensions"][[3]],
		    True,
		    Message[ExportBinaryIBPMatrix::nodims]; Close[stream]; Return[$Failed]
		];

        (* incre \:53ef\:80fd\:76f4\:63a5\:5b58\:50a8\:5728 regData \:4e2d\:ff0c\:82e5\:65e0\:5219\:5c1d\:8bd5\:4ece\:67d0\:4e2a\:7b97\:5b50\:63a8\:65ad\:ff1f\:9ed8\:8ba4\:4e3a 2 *)
        incre = Lookup[regData, "incre", 2];

        (* \:5199\:5165\:5143\:6570\:636e *)
        BinaryWrite[stream, {nibp, ne, nb, incre}, "Integer32", ByteOrdering -> +1];
        BinaryWrite[stream, mod, "Integer64", ByteOrdering -> +1];

        (* \:5199\:5165\:6bcf\:4e2a\:7b97\:5b50 *)
        Scan[
            Function[opName,
                exists = KeyExistsQ[regData, opName];
                BinaryWrite[stream, Boole[exists], "Byte", ByteOrdering -> +1];
                If[exists,
                    opData = regData[opName];
                    If[Head[opData] === SparseArray,
                        binaryData =sparseToBinary[opData];
                       {dimsLen, dims, rowPtrLen, rowPtr, colIdxLen, colIdx, valsLen, vals} = binaryData;
                        typeList = Join[
						    ConstantArray["Integer32", 1 + dimsLen + 1 + rowPtrLen + 1 + colIdxLen + 1],
						    If[mod == 0, ConstantArray["Real64", valsLen], ConstantArray["Integer64", valsLen]]
						];
						BinaryWrite[stream, binaryData, typeList, ByteOrdering -> +1]
						,(*\:5982\:679c\:4e0d\:662f SparseArray\:ff0c\:6253\:5370\:8b66\:544a\:5e76\:8df3\:8fc7 *)
                        Print["Warning: ", opName, " is not a SparseArray, skipping."];
                        (* \:4f46\:4e3a\:4e86\:4fdd\:6301\:683c\:5f0f\:4e00\:81f4\:ff0c\:4ecd\:9700\:5199\:5165\:7a7a\:6570\:636e\:ff1f\:5e94\:8be5\:8df3\:8fc7\:ff0c\:4f46\:4e4b\:524d\:5df2\:7ecf\:5199\:5165 exists=1\:ff0c\:6240\:4ee5\:5fc5\:987b\:5199\:5165\:6570\:636e\:3002\:6211\:4eec\:5c1d\:8bd5\:5c06\:5176\:8f6c\:6362\:4e3a\:7a00\:758f\:ff1f\:6216\:5199\:5165\:7a7a\:7a00\:758f\:3002 *)
                        (* \:6b64\:5904\:7b80\:5316\:ff1a\:5982\:679c\:975e\:7a00\:758f\:ff0c\:6211\:4eec\:5199\:4e00\:4e2a\:7a7a\:7684\:7a00\:758f\:7ed3\:6784\:ff08dims_len=0\:ff09*)
                        BinaryWrite[stream, {0}, "Integer32", ByteOrdering -> +1]; (* dims_len=0 *)
                    ]
                ]
            ],
            operatorOrder
        ],
        {i, numRegs}
    ];

    Close[stream];
    Print["Binary export completed: ", filename];
    filename
];

(* \:8f85\:52a9\:51fd\:6570\:ff1a\:8ba1\:7b97\:5355\:4e2a Ring \:7684\:6240\:6709\:6570\:636e *)
ComputeRingData[ringData_, Alist_, ne_, modulus_] := Module[
    {ringdef, ringmat, varRule, varIndep, varIndepSym, varConvRule, basisIdx, coeffs, matrices, matricesInv, limitSector, dims, flatA, flatAinv},
    (* 1. \:63d0\:53d6\:57fa\:7840\:5b9a\:4e49 *)
    ringdef = ringData["CoordinateRing"];
    limitSector = ringdef["LimitSector"]; (* Integer List *)
    (* 2. \:6784\:5efa\:57fa\:77e9\:9635 (M) *)
    ringmat = N@Values[ringdef["MonomialBasisMatrix"]];
    (* 3. \:9884\:63d0\:53d6\:5c5e\:6027 *)
    varRule = ringdef["VarRule"];
    varIndep = ringdef["VarIndep"];
    (* \:5982\:679c varIndep \:662f\:5b57\:7b26\:4e32\:5217\:8868\:ff0c\:8f6c\:6362\:4e3a\:7b26\:53f7 *)
    If[Length[varIndep] > 0 && StringQ[varIndep[[1]]],
        varIndepSym = Table[Unique["a"], {Length[varIndep]}];
        varConvRule = Thread[varIndep -> varIndepSym];
        varIndep = varIndepSym,
        varConvRule = {}
    ];
    basisIdx = ringdef["MonomialBasisIndex"];
    (* 4. \:8ba1\:7b97\:7cfb\:6570\:77e9\:9635 (C) *)
    coeffs = Map[If[varIndep =!= {}, Lookup[monomialRulesPower[# /. varRule /. varConvRule, varIndep], basisIdx, 0], {# /. varRule}]&, Alist[[1 ;; ne]]];
    (* 5. \:751f\:6210\:77e9\:9635\:5217\:8868 A_list = C.M *)
    matrices = PolynomialMod[coeffs . ringmat, modulus];
    (* 6. \:8ba1\:7b97\:9006\:77e9\:9635 Ainv_list *)
    matricesInv = Map[Inverse[#,Modulus->modulus]&, matrices];
    (* 7. \:6570\:636e\:6574\:5f62\:4e3a C++ \:53cb\:597d\:7684\:683c\:5f0f *)
    flatA = Flatten[#, 1] & /@ matrices;
    flatAinv = Flatten[#, 1] & /@ matricesInv;
    (* \:8fd4\:56de\:6570\:636e\:5305 *)
    <|"LimitSector" -> limitSector, "FlatA" -> flatA, "FlatAinv" -> flatAinv, "MatDim" -> Length[ringmat]|>
];

(* \:4e3b\:5bfc\:51fa\:51fd\:6570\:ff1aRing \:6570\:636e *)
ExportBinaryRingData[filename_String, expregdata_, Alist_, ne_, modulus_] := Module[
    {stream, flatExpreg, count, cellData, matCount, matSize, nb, dataType},
    flatExpreg = Flatten[Values@expregdata];
    count = Length[flatExpreg];
    stream = OpenWrite[filename, BinaryFormat -> True];
    If[stream === $Failed, Return[$Failed]];

    (* \:6587\:4ef6\:5934 *)
    BinaryWrite[stream, count, "Integer32"];  (* \:5199\:5165\:603b\:7684\:6570\:636e\:5355\:5143(Ring)\:6570\:91cf *)
    BinaryWrite[stream, ne, "Integer32"];     (* \:5199\:5165\:6bcf\:4e2a\:5355\:5143\:5305\:542b\:7684\:77e9\:9635\:6570\:91cf *)
    Print["Exporting ", count, " cells to ", filename, "..."];

    Do[
        cellData = ComputeRingData[flatExpreg[[i]], Alist, ne, modulus];
        (* 1. \:5199\:5165 LimitSector\:ff1a\:5148\:5199\:957f\:5ea6\:ff0c\:518d\:5199\:6570\:636e *)
        BinaryWrite[stream, Length[cellData["LimitSector"]], "Integer32"];
        BinaryWrite[stream, cellData["LimitSector"], "Integer32"];
        (* 2. \:5199\:5165\:77e9\:9635\:7ef4\:5ea6\:4fe1\:606f (nb) *)
        nb = cellData["MatDim"];
        BinaryWrite[stream, nb, "Integer32"];
        If[modulus==0, dataType="Real64", dataType="Integer64"];
        (* 3. \:5199\:5165 A_list (Flattened) *)
        BinaryWrite[stream, Flatten[cellData["FlatA"]], dataType];
        (* 4. \:5199\:5165 Ainv_list (Flattened) *)
        BinaryWrite[stream, Flatten[cellData["FlatAinv"]], dataType];
        , {i, count}
    ];

    Close[stream];
    Print["Export Complete: ", filename];
    filename
];

End[]
EndPackage[]
