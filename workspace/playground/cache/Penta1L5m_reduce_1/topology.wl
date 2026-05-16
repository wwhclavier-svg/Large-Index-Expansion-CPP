If[KeyExistsQ["ParallelOptions"/.SystemOptions["ParallelOptions"],#],SetSystemOptions["ParallelOptions"->#->1]]&/@{"ParallelThreadNumber","MKLThreadNumber"};
SetDim[4 - 2*eps];
Internal = {l1};
External = {p1, p2, p3, p4};
MomentumConservation = {};
Replacements = {p1^2 -> 0, p1*p2 -> s12/2, p1*p3 -> (-s12 - s23 + s45)/2, p1*p4 -> (-s15 + s23 - s45)/2, p2^2 -> 0, p2*p3 -> s23/2, p2*p4 -> (s15 - s23 - s34)/2, p3^2 -> 0, p3*p4 -> s34/2, p4^2 -> 0};
Propagators = {-l1^2 + msq, msq - (l1 + p1)^2, msq - (l1 + p1 + p2)^2, msq - (l1 + p1 + p2 + p3)^2, msq - (l1 + p1 + p2 + p3 + p4)^2}/. MomentumConservation//Expand;
NewBasis[Penta1L5m,Propagators,Internal, External, MomentumConservation, Replacements,"ExtraIntDeriv"->{{{0}, {0}, {0}, {0}, {0}}}];
AnalyzeSectors[Penta1L5m,{1, 1, 1, 1, 1},"CutDs"->{0, 0, 0, 0, 0},"Prescription"->{1, 1, 1, 1, 1},"CloseSyms"->False,"ExtSyms"->True];
