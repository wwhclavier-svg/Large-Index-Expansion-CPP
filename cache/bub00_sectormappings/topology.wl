If[KeyExistsQ["ParallelOptions"/.SystemOptions["ParallelOptions"],#],SetSystemOptions["ParallelOptions"->#->1]]&/@{"ParallelThreadNumber","MKLThreadNumber"};
SetDim[4 - 2*eps];
Internal = {k1};
External = {p1};
MomentumConservation = {};
Replacements = {p1^2 -> s};
Propagators = {-k1^2, -(k1 - p1)^2}/. MomentumConservation//Expand;
NewBasis[bub00,Propagators,Internal, External, MomentumConservation, Replacements,"ExtraIntDeriv"->{{{0}, {0}}}];
AnalyzeSectors[bub00,{1, 1},"CutDs"->{0, 0},"Prescription"->{1, 1},"CloseSyms"->False,"ExtSyms"->True];
