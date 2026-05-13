If[KeyExistsQ["ParallelOptions"/.SystemOptions["ParallelOptions"],#],SetSystemOptions["ParallelOptions"->#->1]]&/@{"ParallelThreadNumber","MKLThreadNumber"};
SetDim[4 - 2*eps];
Internal = {l1, l2};
External = {p};
MomentumConservation = {};
Replacements = {p^2 -> s};
Propagators = {-l1^2 - msq, -(l1 + p)^2, -l2^2 - msq, -(l2 + p)^2, -msq - (l1 + l2 + p)^2}/. MomentumConservation//Expand;
NewBasis[SR212m3,Propagators,Internal, External, MomentumConservation, Replacements,"ExtraIntDeriv"->{{{0}, {0}, {0}, {0}, {0}}}];
AnalyzeSectors[SR212m3,{1, 0, 1, 0, 1},"CutDs"->{0, 0, 0, 0, 0},"Prescription"->{1, 1, 1, 1, 1},"CloseSyms"->False,"ExtSyms"->True];
