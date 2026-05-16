(* ::Package:: *)

BeginPackage["ImportIBPMatrixBinary`"]

ImportBinaryIBPMatrix::usage = "ImportBinaryIBPMatrix[filename] reads an IBPMat_*.bin file (big-endian) and returns the region data."
ImportBinaryRingData::usage = "ImportBinaryRingData[filename] reads a RingData_*.bin file (little-endian) and returns the cell data."

Begin["`Private`"]

(* Helper: read a single value from stream with given format and endianness *)
readOne[stream_, fmt_, endian_] := BinaryReadList[stream, fmt, 1, ByteOrdering -> endian][[1]];

(* ================================================================ *)
(* ImportBinaryIBPMatrix
   Format: big-endian (written by MMA ExportIBPMatrixBinary` with ByteOrdering->+1)
   Magic "IBP1" → numRegs(int32) → per region:
     nibp,ne,nb,incre(int32) mod(int64)
     8 operator flags (byte)
       per existing op: dimsLen(int32), dims[int32], rowPtrLen(int32), rowPtr[int32],
                        colIdxLen(int32), colIdx[int32], valsLen(int32), vals[int64]
*)
(* ================================================================ *)
ImportBinaryIBPMatrix[filename_String] := Module[
    {stream, magic, numRegs, regions, nibp, ne, nb, incre, mod,
     operatorOrder, exists, dimsLen, dims, rowPtrLen, rowPtr, rowIndices,
     colIdxLen, colIdx, valsLen, vals, regData, i, dataType,
     ndims, colPositions, positions},

    stream = OpenRead[filename, BinaryFormat -> True];
    If[stream === $Failed, Print["ERROR: Cannot open ", filename]; Return[$Failed]];

    magic = BinaryReadList[stream, "Character8", 4];
    If[magic =!= {"I", "B", "P", "1"},
        Print["ERROR: Invalid magic number: ", StringJoin[magic]];
        Close[stream]; Return[$Failed]
    ];
    numRegs = readOne[stream, "Integer32", +1];
    Print["Reading ", numRegs, " regions from ", filename];

    operatorOrder = {"M1", "N1", "K1", "F0", "F2", "K1s", "K2s", "F2s"};
    regions = Table[
        nibp  = readOne[stream, "Integer32", +1];
        ne    = readOne[stream, "Integer32", +1];
        nb    = readOne[stream, "Integer32", +1];
        incre = readOne[stream, "Integer32", +1];
        mod   = readOne[stream, "Integer64", +1];
        dataType = If[mod == 0, "Real64", "Integer64"];

        regData = <|"nibp" -> nibp, "ne" -> ne, "nb" -> nb, "incre" -> incre, "mod" -> mod|>;

        Scan[
            Function[opName,
                exists = readOne[stream, "Byte", +1];
                If[exists === 1,
                    dimsLen = readOne[stream, "Integer32", +1];
                    If[dimsLen > 0,
                        dims = BinaryReadList[stream, "Integer32", dimsLen, ByteOrdering -> +1];
                        rowPtrLen = readOne[stream, "Integer32", +1];
                        rowPtr = BinaryReadList[stream, "Integer32", rowPtrLen, ByteOrdering -> +1];
                        colIdxLen = readOne[stream, "Integer32", +1];
                        colIdx = BinaryReadList[stream, "Integer32", colIdxLen, ByteOrdering -> +1];
                        valsLen = readOne[stream, "Integer32", +1];
                        vals = BinaryReadList[stream, dataType, valsLen, ByteOrdering -> +1];
                        rowIndices = Join @@ MapThread[ConstantArray, {Range[Length[rowPtr] - 1], Differences[rowPtr]}];
                        ndims = Length[dims];
                        If[ndims == 2,
                            regData[opName] = SparseArray[Transpose[{rowIndices, colIdx}] -> vals, dims],
                            colPositions = Partition[colIdx, ndims - 1];
                            positions = Transpose[Prepend[Transpose[colPositions], rowIndices]];
                            regData[opName] = SparseArray[positions -> vals, dims]
                        ];
                    ];
                ];
            ],
            operatorOrder
        ];
        regData
        , {i, numRegs}
    ];

    Close[stream];
    Print["Done: ", numRegs, " regions loaded."];
    regions
];

(* ================================================================ *)
(* ImportBinaryRingData
   Format: little-endian (written by C++ BinaryRingWriter.hpp, native endian)
   Header: count(int32) ne(int32)
   Per cell: secLen(int32) limitSector[int32×secLen] nb(int32)
             A_flat[int64×ne×nb×nb] Ainv_flat[int64×ne×nb×nb]
*)
(* ================================================================ *)
ImportBinaryRingData[filename_String] := Module[
    {stream, count, ne, cells, lsLen, limitSector, nb, flatA, flatAinv,
     totalA, dataType, i},

    stream = OpenRead[filename, BinaryFormat -> True];
    If[stream === $Failed, Print["ERROR: Cannot open ", filename]; Return[$Failed]];

    count = readOne[stream, "Integer32", -1];
    ne    = readOne[stream, "Integer32", -1];
    Print["Reading ", count, " cells, ne=", ne, " from ", filename];

    dataType = "Integer64";
    cells = Table[
        lsLen = readOne[stream, "Integer32", -1];
        limitSector = BinaryReadList[stream, "Integer32", lsLen, ByteOrdering -> -1];
        nb = readOne[stream, "Integer32", -1];
        totalA = ne * nb * nb;

        flatA = BinaryReadList[stream, dataType, totalA, ByteOrdering -> -1];
        flatAinv = BinaryReadList[stream, dataType, totalA, ByteOrdering -> -1];

        <|
            "LimitSector" -> limitSector,
            "nb" -> nb,
            "FlatA" -> Table[Partition[flatA[[(i-1)*nb*nb + 1 ;; i*nb*nb]], nb], {i, ne}],
            "FlatAinv" -> Table[Partition[flatAinv[[(i-1)*nb*nb + 1 ;; i*nb*nb]], nb], {i, ne}]
        |>
        , {i, count}
    ];

    Close[stream];
    Print["Done: ", count, " cells loaded."];
    cells
];

End[]
EndPackage[]
