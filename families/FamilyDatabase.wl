(* ::Package:: *)

(* ============================================================ *)
(* FamilyDatabase.wl \[LongDash] IBP \:79ef\:5206\:65cf\:96c6\:4e2d\:5b9a\:4e49\:6570\:636e\:5e93 (+ JSON \:540c\:6b65\:5668)    *)
(* ============================================================ *)
(* \:7528\:6cd5:                                                           *)
(*   Get["families/FamilyDatabase.wl"];               *)
(*   fam = "bub00";                                                *)
(*   config = $FamilyDatabase[fam];                                *)
(*   LIEDefineFamily[                                              *)
(*       config["Propagators"], config["LoopMomenta"],             *)
(*       config["ExternalMomenta"], config["KinematicRules"],      *)
(*       config["TopSector"], "Numeric" -> config["Numeric"],      *)
(*       Modulus -> config["Modulus"]]                             *)
(* ============================================================ *)

$FamilyDatabase = <|

    (* ======================================================== *)
    (* 1-Loop \:65cf (L=1)                                          *)
    (* ======================================================== *)

    (* -------------------------------------------------------- *)
    (* bub00 \[LongDash] 1-loop bubble, massless (msq=0)                  *)
    (* \:4f20\:64ad\:5b50: {-k1^2, -(k1-p1)^2}                               *)
    (* -------------------------------------------------------- *)
    "bub00" -> <|
        "Description" -> "1-loop bubble, msq=0",
        "Propagators" -> {k1^2, (k1 - p1)^2},
        "LoopMomenta" -> {k1},
        "ExternalMomenta" -> {p1},
        "KinematicRules" -> {p1^2 -> s},
        "TopSector" -> {1, 1},
        "Numeric" -> {s -> 3, msq -> 0, "d" -> 1/3},
        "Modulus" -> Prime[10000000],  (* 179424673 *)
        "Graph"->{{1,2},{2,1}}
    |>,

    (* -------------------------------------------------------- *)
    (* bub10 \[LongDash] 1-loop bubble, 1-massive (msq=1)                 *)
    (* \:4f20\:64ad\:5b50: {-k1^2+msq, -(k1-p1)^2}                           *)
    (* -------------------------------------------------------- *)
    "bub10" -> <|
        "Description" -> "1-loop bubble, msq=1",
        "Propagators" -> {k1^2 - msq, (k1 - p1)^2},
        "LoopMomenta" -> {k1},
        "ExternalMomenta" -> {p1},
        "KinematicRules" -> {p1^2 -> s},
        "TopSector" -> {1, 1},
        "Numeric" -> {s -> 3, msq -> 1, "d" -> 1/3},
        "Modulus" -> Prime[10000000],
        "Graph"->{{1,2},{2,1}}
    |>,

    (* -------------------------------------------------------- *)
    (* bub11 \[LongDash] 1-loop bubble, both-massive (msq=1)              *)
    (* \:4f20\:64ad\:5b50: {-k1^2 + msq, -(k1-p1)^2 + msq}                   *)
    (* -------------------------------------------------------- *)
    "bub11" -> <|
        "Description" -> "1-loop bubble, msq=1",
        "Propagators" -> {k1^2 - msq, (k1 - p1)^2 - msq},
        "LoopMomenta" -> {k1},
        "ExternalMomenta" -> {p1},
        "KinematicRules" -> {p1^2 -> s},
        "TopSector" -> {1, 1},
        "Numeric" -> {s -> 3, msq -> 1, "d" -> 1/3},
        "Modulus" -> Prime[10000000],
        "Graph"->{{1,2},{2,1}}
    |>,

    (* -------------------------------------------------------- *)
    (* Tri \[LongDash] 1-loop triangle (3 propagators, 2 ext momenta)     *)
    (* Kinematics: p1^2=s1, p2^2=s2, (p1+p2)^2=s3               *)
    (* -------------------------------------------------------- *)
    "Tri" -> <|
        "Description" -> "1-loop triangle, 3 props, s1=3,s2=8,s3=0",
        "Propagators" -> {k1^2 - msq, (k1 - p1)^2 - msq, (k1 + p2)^2 - msq},
        "LoopMomenta" -> {k1},
        "ExternalMomenta" -> {p1, p2},
        "KinematicRules" -> {
            p1^2 -> s1, p2^2 -> s2,
            (p1 + p2)^2 -> s3,
            p1*p2 -> (s3 - s1 - s2)/2
        },
        "TopSector" -> {1, 1, 1},
        "Numeric" -> {s1 -> 3, s2 -> 8, s3 -> 0, msq -> 0, "d" -> 1/3},
        "Modulus" -> Prime[10000000],
        "Graph"->{{1,2},{3,1},{2,3}}
    |>,

    (* -------------------------------------------------------- *)
    (* Box \[LongDash] 1-loop box (4 propagators, 3 ext momenta)          *)
    (* Kinematics: p1^2=s, p2^2=t, p3^2=u, Mandelstam          *)
    (* -------------------------------------------------------- *)
    "Box" -> <|
        "Description" -> "1-loop box, 4 props, s=3,t=2,u=-5",
        "Propagators" -> {
            k1^2 - msq,
            (k1 - p1)^2 - msq,
            (k1 - p1 - p2)^2 - msq,
            (k1 - p1 - p2 - p3)^2 - msq
        },
        "LoopMomenta" -> {k1},
        "ExternalMomenta" -> {p1, p2, p3},
        "KinematicRules" -> {
            p1^2 -> s, p2^2 -> t, p3^2 -> u,
            (p1 + p2)^2 -> 0, (p2 + p3)^2 -> 0,
            (p1 + p2 + p3)^2 -> s + t + u,
            p1*p2 -> (- s - t)/2,
            p2*p3 -> (- t - u)/2,
            p1*p3 ->  (2 t - s - u)/2
        },
        "TopSector" -> {1, 1, 1, 1},
        "Numeric" -> {s -> 3, t -> 2, msq -> 0, "d" -> 1/3, u -> -5},
        "Modulus" -> Prime[10000000],
        "Graph"->{{1,4},{2,1},{3,2},{4,3}}
    |>,

    (* ======================================================== *)
    (* 2-Loop \:65cf (L=2)                                          *)
    (* ======================================================== *)

    (* -------------------------------------------------------- *)
    (* SR \[LongDash] Sunrise (2L1P, 5 propagators)                       *)
    (* \:6765\:6e90: LIETest_2L2P_SR212.nb "Massless" \:8282                *)
    (* -------------------------------------------------------- *)
    "SR212" -> <|
        "Description" -> "Sunrise 2L1P, massless, s=1",
        "Propagators" -> {
            l1^2,
            (l1 + p)^2,
            l2^2,
            (l2 + p)^2,
            (l1 + l2 + p)^2
        },
        "LoopMomenta" -> {l1, l2},
        "ExternalMomenta" -> {p},
        "KinematicRules" -> {p^2 -> s},
        "TopSector" -> {1, 1, 1, 1, 1},
        "Numeric" -> {s -> 1, "d" -> 1/13},
        "Modulus" -> Prime[10000000],
        "Graph"->{{3,1},{1,4},{2,4},{3,2},{4,3}}
    |>,

    (* SR212-3m \[LongDash] Sunrise 3-Massive (2L1P, 5 prop)                  *)
    (* \:6765\:6e90: LIETest_2L2P_SR212_3Massive.nb                     *)
    (* -------------------------------------------------------- *)
    "SR212-3m" -> <|
        "Description" -> "Sunrise 2L1P, 3-massive (props 1,3,5), s=3/2",
        "Propagators" -> {
            l1^2 - msq,
            (l1 + p)^2,
            l2^2 - msq,
            (l2 + p)^2,
            (l1 + l2 + p)^2 - msq
        },
        "LoopMomenta" -> {l1, l2},
        "ExternalMomenta" -> {p},
        "KinematicRules" -> {p^2 -> s},
        "TopSector" -> {1, 0, 1, 0, 1},
        "Numeric" -> {s -> 3/2, msq -> 1, "d" -> 1/7},
        "Modulus" -> Prime[10000000],
        "Graph"->{{3,1},{1,4},{2,4},{3,2},{4,3}}
    |>,

    (* -------------------------------------------------------- *)
    (* SR212-5m \[LongDash] Sunrise All-Massive (2L1P, 5 prop)                *)
    (* \:6765\:6e90: LIETest_2L2P_SR212.nb "All Massive" \:8282             *)
    (* -------------------------------------------------------- *)
    "SR212-5m" -> <|
        "Description" -> "Sunrise 2L1P, all-massive, s=0",
        "Propagators" -> {
            l1^2 - msq,
            (l1 + p)^2 - msq,
            l2^2 - msq,
            (l2 + p)^2 - msq,
            (l1 + l2 + p)^2 - msq
        },
        "LoopMomenta" -> {l1, l2},
        "ExternalMomenta" -> {p},
        "KinematicRules" -> {p^2 -> s},
        "TopSector" -> {1, 1, 1, 1, 1},
        "Numeric" -> {s -> 0, msq -> 1, "d" -> 1/13},
        "Modulus" -> Prime[10000000],
        "Graph"->{{3,1},{1,4},{2,4},{3,2},{4,3}}
    |>,

    (* -------------------------------------------------------- *)
    (* NP222 \[LongDash] Non-Planar 2-loop 2-point (2L2P, 7 propagators)  *)
    (* \:6765\:6e90: LIETest_2L2P_SR212.nb "Non Planar 2-loop" \:8282      *)
    (* -------------------------------------------------------- *)
    "NP222" -> <|
        "Description" -> "Non-Planar 2L2P, s=1,s2=0,m=0",
        "Propagators" -> {
            (l2 + p1)^2 + msq,
            (l1 - l2 - p1 - p2)^2,
            l2^2 + msq,
            (l1 - p2)^2 + msq,
            (l1 - l2)^2,
            l1^2 + msq,
            (l1 + p1)^2
        },
        "LoopMomenta" -> {l1, l2},
        "ExternalMomenta" -> {p1, p2},
        "KinematicRules" -> {
            p1^2 -> 0, p2^2 -> s2,
            p1*p2 -> (s - s2)/2
        },
        "TopSector" -> {1, 1, 1, 1, 1, 1, 0},
        "Numeric" -> {s -> 1, s2 -> 0, msq -> 0, "d" -> 1/13},
        "Modulus" -> Prime[10000000],
        "Graph"->{{1,4},{3,4},{5,1},{4,2},{5,3},{2,5}}
    |>,

    (* -------------------------------------------------------- *)
    (* NP222var \[LongDash] arXiv:2305.13951 non-planar triangle (2L2P)    *)
    (* D1=(l1-p1)^2, D2=(l2-p1)^2-m^2, D3=(l1+p2)^2,            *)
    (* D4=(l1-l2+p2)^2-m^2, D5=(l1-l2)^2-\[Kappa]m^2, D6=l2^2-\[Kappa]m^2,    *)
    (* D7=l1^2 (ISP). \[Kappa]=1\[RightArrow]fam(a) 11MIs, \[Kappa]=0\[RightArrow]fam(b) 18MIs.       *)
    (* \:6765\:6e90: arXiv:2305.13951 Eq.(2), Jiang/Wang/Yang/Zhao 2023  *)
    (* -------------------------------------------------------- *)
    "NP222var" -> <|
        "Description" -> "NP 2L2P triangle, arXiv:2305.13951, k=0(fam b)",
        "Propagators" -> {
            (l1 - p1)^2,
            (l2 - p1)^2 - msq,
            (l1 + p2)^2,
            (l1 - l2 + p2)^2 - msq,
            (l1 - l2)^2 - k*msq,
            l2^2 - k*msq,
            l1^2
        },
        "LoopMomenta" -> {l1, l2},
        "ExternalMomenta" -> {p1, p2},
        "KinematicRules" -> {
            p1^2 -> 0, p2^2 -> 0,
            p1*p2 -> s/2
        },
        "TopSector" -> {1, 1, 1, 1, 1, 1, 0},
        "Numeric" -> {s -> 1, msq -> 0, k -> 0, "d" -> 1/13},
        "Modulus" -> Prime[10000000],
        "Graph"->{{4,3},{4,1},{5,3},{2,5},{4,2},{1,5}}
    |>,

    (* -------------------------------------------------------- *)
    (* TB123 \[LongDash] Tri-Box / Tennis Ball 2L3P (7 prop, 6\:4e3b+1\:8f85)     *)
    (* \:6765\:6e90: LIETest_2L3P_TB123&NP222_Eigenvalue.nb Cell 3,5   *)
    (* -------------------------------------------------------- *)
    "TB123" -> <|
        "Description" -> "Tri-Box 2L3P, massless, s=3,s1=2,s2=4",
        "Propagators" -> {
            l1^2,                         (* 1 *)
            (l1 - p1)^2,                  (* 2 *)
            (l1 - p1 - p2)^2,             (* 3 *)
            (l2 + p1 + p2)^2,             (* 4 *)
            l2^2,                         (* 5 *)
            (l1 + l2)^2,                  (* 6 *)
            (l2 + p1)^2                   (* 7: \:8f85\:52a9 ISP *)
        },
        "LoopMomenta" -> {l1, l2},
        "ExternalMomenta" -> {p1, p2},
        "KinematicRules" -> {
            p1^2 -> s1, p2^2 -> s2,
            p1*p2 -> (s - s1 - s2)/2
        },
        "TopSector" -> {1, 1, 1, 1, 1, 1, 0},
        "Numeric" -> {s -> 3, s1 -> 2, s2 -> 4, msq -> 0, "d" -> 1/13},
        "Modulus" -> Prime[10000000],
        "Graph"->{{1,4},{2,1},{5,2},{5,3},{3,4},{4,5}}
    |>,

    (* -------------------------------------------------------- *)
    (* TB123m \[LongDash] Tri-Box / Tennis Ball 2L3P massive              *)
    (* \:6765\:6e90: LIETest_2L3P_TB123&NP222_Eigenvalue.nb Cell 9,11  *)
    (* -------------------------------------------------------- *)
    "TB123m" -> <|
        "Description" -> "Tri-Box 2L3P, massive, s=3,s1=2,s2=5,m=1",
        "Propagators" -> {
            l1^2 - msq,
            (l1 - p1)^2 - msq,
            (l1 - p1 - p2)^2 - msq,
            (l2 + p1 + p2)^2 - msq,
            l2^2 - msq,
            (l1 + l2)^2 - msq,
            (l2 + p1)^2
        },
        "LoopMomenta" -> {l1, l2},
        "ExternalMomenta" -> {p1, p2},
        "KinematicRules" -> {
            p1^2 -> s1, p2^2 -> s2,
            p1*p2 -> (s - s1 - s2)/2
        },
        "TopSector" -> {1, 1, 1, 1, 1, 1, 0},
        "Numeric" -> {s -> 3, s1 -> 2, s2 -> 5, msq -> 1, "d" -> 1/13},
        "Modulus" -> Prime[10000000],
        "Graph"->{{1,4},{2,1},{5,2},{5,3},{3,4},{4,5}}
    |>,

    (* -------------------------------------------------------- *)
    (* DB313 \[LongDash] Double Box (2L4P, 9 prop, 7\:4e3b+2\:8f85)              *)
    (* \:6765\:6e90: LIETest_2L4P_DB313_New.nb, dbox_data.wl           *)
    (* Kinematics: s=(p1+p2)^2, t=(p1+p4)^2, massless ext       *)
    (* -------------------------------------------------------- *)
    "DB313" -> <|
        "Description" -> "Double Box 2L4P, s=2,t=3",
        "Propagators" -> {
            l1^2,                   (* 1 *)
            (l1 + p1)^2,            (* 2 *)
            (l1 + p1 + p2)^2,       (* 3 *)
            (l2 - p1 - p2)^2,       (* 4 *)
            (l2 + p4)^2,            (* 5 *)
            l2^2,                   (* 6 *)
            (l1 + l2)^2,            (* 7 *)
            (l2 + p2)^2,            (* 8: \:8f85\:52a9 *)
            (l1 + p4)^2             (* 9: \:8f85\:52a9 *)
        },
        "LoopMomenta" -> {l1, l2},
        "ExternalMomenta" -> {p1, p2, p4},
        "KinematicRules" -> {
            p1^2 -> 0, p2^2 -> 0, p4^2 -> 0,
            p1*p2 -> s/2,
            p2*p4 -> -(s + t)/2,
            p1*p4 -> t/2
        },
        "TopSector" -> {1, 1, 1, 1, 1, 1, 1, 0, 0},
        "Numeric" -> {s -> 2, t -> 3, "d" -> 1/3},
        "Modulus" -> Prime[10000000],
        "Graph"->{{5,1},{1,2},{2,6},{3,6},{4,3},{5,4},{6,5}}
    |>,

    (* -------------------------------------------------------- *)
    (* DB313-4m \[LongDash] Double Box 2L4P, 4-massive (props 1,2,3,7)     *)
    (* \:57fa\:4e8e DB313\:ff0c\:5728 propagators 0,1,2,6 (1-indexed: 1,2,3,7)   *)
    (* \:4e0a\:52a0 msq.                                                   *)
    (* -------------------------------------------------------- *)
    "DB313-4m" -> <|
        "Description" -> "Double Box 2L4P, 4-massive (props 1,2,3,7), s=2,t=3,msq=1",
        "Propagators" -> {
            l1^2 + msq,             (* 1 +msq *)
            (l1 + p1)^2 + msq,      (* 2 +msq *)
            (l1 + p1 + p2)^2 + msq, (* 3 +msq *)
            (l2 - p1 - p2)^2,       (* 4 *)
            (l2 + p4)^2,            (* 5 *)
            l2^2,                   (* 6 *)
            (l1 + l2)^2 + msq,      (* 7 +msq *)
            (l2 + p2)^2,            (* 8: \:8f85\:52a9 *)
            (l1 + p4)^2             (* 9: \:8f85\:52a9 *)
        },
        "LoopMomenta" -> {l1, l2},
        "ExternalMomenta" -> {p1, p2, p4},
        "KinematicRules" -> {
            p1^2 -> 0, p2^2 -> 0, p4^2 -> 0,
            p1*p2 -> s/2,
            p2*p4 -> -(s + t)/2,
            p1*p4 -> t/2
        },
        "TopSector" -> {1, 1, 1, 1, 1, 1, 1, 0, 0},
        "Numeric" -> {s -> 2, t -> 3, msq -> 1, "d" -> 1/3},
        "Modulus" -> Prime[10000000],
        "Graph"->{{5,1},{1,2},{2,6},{3,6},{4,3},{5,4},{6,5}}
    |>,

    (* -------------------------------------------------------- *)
    (* NP322 \[LongDash] Non-Planar 2L4P (9 propagators)                  *)
    (* \:6765\:6e90: LIETest_2L4P_NP322.nb Cell 1-3                     *)
    (* \:4e0e DB313 (Double Box) \:4e0d\:540c: NP322 \:662f\:975e\:5e73\:9762\:5bf9\:5e94\:7269          *)
    (* -------------------------------------------------------- *)
    "NP322" -> <|
        "Description" -> "Non-Planar 2L4P, massless, s=1,t=5",
        "Propagators" -> {
            (k2 + p1)^2 - msq,             (* 1 *)
            (k1 - k2 - p1 - p2)^2,         (* 2 *)
            k2^2 - msq,                    (* 3 *)
            (k1 - p2)^2 - msq,             (* 4 *)
            (k1 - k2)^2,                   (* 5 *)
            k1^2 - msq,                    (* 6 *)
            (k1 - k2 + p4)^2,              (* 7 *)
            (k1 + p1)^2,                   (* 8 *)
            (k1 + p4)^2                    (* 9 *)
        },
        "LoopMomenta" -> {k1, k2},
        "ExternalMomenta" -> {p1, p2, p4},
        "KinematicRules" -> {
            p1^2 -> m1, p2^2 -> m2, p4^2 -> 0,
            p1*p2 -> (s - m1 - m2)/2,
            p1*p4 -> (t - m1)/2,
            p2*p4 -> -(s + t - m1)/2
        },
        "TopSector" -> {1, 1, 1, 1, 1, 1, 1, 0, 0},
        "Numeric" -> {s -> 1, t -> 5, m1 -> 0, m2 -> 0, msq -> 0, "d" -> 1/7},
        "Modulus" -> Prime[10000000],
        "Graph"->{{1,5},{3,5},{6,1},{5,2},{6,4},{2,6},{4,3}}
    |>,

    (* -------------------------------------------------------- *)
    (* NP322m \[LongDash] Non-Planar 2L4P massive variant                  *)
    (* \:6765\:6e90: LIETest_2L4P_NP322.nb Cell 1 (massive numeric)     *)
    (* -------------------------------------------------------- *)
    "NP322m" -> <|
        "Description" -> "Non-Planar 2L4P, massive, s=3,t=2,m1=1,m2=4,m=1",
        "Propagators" -> {
            (k2 + p1)^2 - msq,
            (k1 - k2 - p1 - p2)^2,
            k2^2 - msq,
            (k1 - p2)^2 - msq,
            (k1 - k2)^2,
            k1^2 - msq,
            (k1 - k2 + p4)^2,
            (k1 + p1)^2,
            (k1 + p4)^2
        },
        "LoopMomenta" -> {k1, k2},
        "ExternalMomenta" -> {p1, p2, p4},
        "KinematicRules" -> {
            p1^2 -> m1, p2^2 -> m2, p4^2 -> 0,
            p1*p2 -> (s - m1 - m2)/2,
            p1*p4 -> (t - m1)/2,
            p2*p4 -> -(s + t - m1)/2
        },
        "TopSector" -> {1, 1, 1, 1, 1, 1, 1, 0, 0},
        "Numeric" -> {s -> 3, t -> 2, m1 -> 1, m2 -> 4, msq -> 1, "d" -> 1/3},
        "Modulus" -> Prime[10000000],
        "Graph"->{{1,5},{3,5},{6,1},{5,2},{6,4},{2,6},{4,3}}
    |>,

    (* -------------------------------------------------------- *)
    (* DP323 \[LongDash] Double Pentagon (2L5P, 11 prop, 8\:4e3b+3\:8f85)        *)
    (* \:6765\:6e90: dpen_data.wl, LIETest_2L5P_DP323.nb                *)
    (* -------------------------------------------------------- *)
    "DP323" -> <|
        "Description" -> "Double Pentagon 2L5P, s12=4,s13=7,s14=8,s23=1,s24=2,m=3",
        "Propagators" -> {
            k1^2 - msq,                            (* 1  *)
            (k1 - p1)^2 - msq,                     (* 2  *)
            (k1 - p1 - p2)^2 - msq,                (* 3  *)
            k2^2 + msq,                            (* 4  *)
            (k2 - p1 - p2 - p3)^2 - msq,           (* 5  *)
            (k2 - p1 - p2 - p3 - p4)^2 - msq,      (* 6  *)
            (k1 - k2)^2 - msq,                     (* 7  *)
            (k1 - k2 + p3)^2 - msq,                (* 8  *)
            (k1 - p1 - p2 - p3 - p4)^2,            (* 9: \:8f85\:52a9 ISP *)
            (k2 - p1)^2,                           (* 10: \:8f85\:52a9 ISP *)
            (k2 - p1 - p2)^2                       (* 11: \:8f85\:52a9 ISP *)
        },
        "LoopMomenta" -> {k1, k2},
        "ExternalMomenta" -> {p1, p2, p3, p4},
        "KinematicRules" -> {
            p1^2 -> 0, p2^2 -> 0, p3^2 -> 0, p4^2 -> 0,
            p1*p2 -> s12/2,
            p1*p3 -> s13/2,
            p1*p4 -> s14/2,
            p2*p3 -> s23/2,
            p2*p4 -> s24/2,
            p3*p4 -> (-s12 - s13 - s14 - s23 - s24)/2
        },
        "TopSector" -> {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0},
        "Numeric" -> {s12 -> 4, s13 -> 7, s14 -> 8, s23 -> 1, s24 -> 2, msq -> 3, "d" -> 1/3},
        "Modulus" -> Prime[10000000],
        "Graph"->{{1,6},{2,1},{7,2},{6,5},{4,7},{5,4},{6,3},{3,7}}
    |>,

    (* ======================================================== *)
    (* 3-Loop \:65cf (L=3)                                          *)
    (* \:6765\:6e90: LIETest_2L4P_TP412&NP322_Eigenvalue.nb             *)
    (* ======================================================== *)

    (* -------------------------------------------------------- *)
    (* BN3L \[LongDash] 3-Loop Banana vacuum (3L0P, 6 propagators)        *)
    (* \:6765\:6e90: TP412&NP322_Eigenvalue.nb Cell 50                   *)
    (* -------------------------------------------------------- *)
    "BN3L" -> <|
        "Description" -> "3-Loop Banana, vacuum, msq=1",
        "Propagators" -> {
            l1^2 - msq,
            l2^2 - msq,
            l3^2 - msq,
            (l1 - l2)^2 - msq,
            (l2 - l3)^2 - msq,
            (l3 - l1)^2 - msq
        },
        "LoopMomenta" -> {l1, l2, l3},
        "ExternalMomenta" -> {},
        "KinematicRules" -> {},
        "TopSector" -> {1, 1, 1, 1, 1, 1},
        "Numeric" -> {msq -> 1, "d" -> 1/13},
        "Modulus" -> Prime[10000000],
        "Graph"->{{1,2},{2,3},{3,1},{2,4},{3,4},{1,4}}
    |>,

    (* -------------------------------------------------------- *)
    (* BN3L1P \[LongDash] 3-Loop Banana 1-Point (9 prop, 8\:4e3b+1\:8f85)        *)
    (* \:6765\:6e90: TP412&NP322_Eigenvalue.nb Cell 56                   *)
    (* -------------------------------------------------------- *)
    "BN3L1P" -> <|
        "Description" -> "3-Loop Banana 1-Point, massless, s=1",
        "Propagators" -> {
            l1^2,                 (* 1 *)
            (l1 + p)^2,           (* 2 *)
            (l2 + p)^2,           (* 3 *)
            l2^2,                 (* 4 *)
            l3^2,                 (* 5 *)
            (l1 - l2)^2,          (* 6 *)
            (l2 - l3)^2,          (* 7 *)
            (l3 - l1)^2,          (* 8 *)
            (l3 - p)^2            (* 9: \:8f85\:52a9 ISP *)
        },
        "LoopMomenta" -> {l1, l2, l3},
        "ExternalMomenta" -> {p},
        "KinematicRules" -> {p^2 -> s},
        "TopSector" -> {1, 1, 1, 1, 1, 1, 1, 1, 0},
        "Numeric" -> {s -> 1, "d" -> 1/13},
        "Modulus" -> Prime[10000000],
        "Graph"->{{3,1},{1,4},{4,2},{2,5},{5,3},{4,6},{5,6},{3,6}}
    |>,

    (* -------------------------------------------------------- *)
    (* BN3L1Pm \[LongDash] 3-Loop Banana 1-Point massive (9 prop, 8\:4e3b+1\:8f85) *)
    (* \:6765\:6e90: TP412&NP322_Eigenvalue.nb Cell 62                   *)
    (* -------------------------------------------------------- *)
    "BN3L1Pm" -> <|
        "Description" -> "3-Loop Banana 1-Point, massive, s=1,msq=1",
        "Propagators" -> {
            l1^2 + msq,
            (l1 + p)^2 + msq,
            (l2 + p)^2 + msq,
            l2^2 + msq,
            l3^2 + msq,
            (l1 - l2)^2 + msq,
            (l2 - l3)^2 + msq,
            (l3 - l1)^2 + msq,
            (l3 - p)^2
        },
        "LoopMomenta" -> {l1, l2, l3},
        "ExternalMomenta" -> {p},
        "KinematicRules" -> {p^2 -> s},
        "TopSector" -> {1, 1, 1, 1, 1, 1, 1, 1, 0},
        "Numeric" -> {s -> 1, msq -> 1, "d" -> 1/13},
        "Modulus" -> Prime[10000000],
        "Graph"->{{3,1},{1,4},{4,2},{2,5},{5,3},{4,6},{5,6},{3,6}}
    |>,

    (* ======================================================== *)
    (* C++ \:8bd5\:9a8c\:6027\:65cf (\:7531 C++ JSON \:5bfc\:5165, 2026-05-13)                *)
    (* ======================================================== *)

    (* -------------------------------------------------------- *)
    (* DP2L5P \[LongDash] Non-Planar DoublePentagon 2L5P, massless        *)
    (* \:6765\:6e90: arXiv:2009.07803 Eq.(3.2) topology (d)              *)
    (* -------------------------------------------------------- *)
    "DP2L5P" -> <|
        "Description" -> "Non-Planar DoublePentagon 2L5P massless, s12=3,s23=2,s34=5,s45=7,s15=11",
        "Propagators" -> {
            l1^2,
            (l1 - p1)^2,
            (l1 - p1 - p2)^2,
            l2^2,
            (l2 - p1 - p2 - p3)^2,
            (l2 - p1 - p2 - p3 - p4)^2,
            (l1 - l2)^2,
            (l1 - l2 + p3)^2,
            (l1 - p1 - p2 - p3 - p4)^2,
            (l2 - p1)^2,
            (l2 - p1 - p2)^2
        },
        "LoopMomenta" -> {l1, l2},
        "ExternalMomenta" -> {p1, p2, p3, p4},
        "KinematicRules" -> {
            p1^2 -> 0, p2^2 -> 0, p3^2 -> 0, p4^2 -> 0,
            p1*p2 -> s12/2, p2*p3 -> s23/2, p3*p4 -> s34/2,
            p1*p3 -> (s45 - s12 - s23)/2,
            p2*p4 -> (s15 - s23 - s34)/2,
            p1*p4 -> (s23 - s45 - s15)/2
        },
        "TopSector" -> {1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0},
        "Numeric" -> {s12 -> 3, s23 -> 2, s34 -> 5, s45 -> 7, s15 -> 11, "d" -> 1/13},
        "Modulus" -> Prime[10000000],
        "Graph"->{{1,6},{2,1},{7,2},{6,5},{4,7},{5,4},{6,3},{3,7}}
    |>,

    (* -------------------------------------------------------- *)
    (* HB2L5P \[LongDash] Non-Planar HexaBox 2L5P, massless               *)
    (* \:6765\:6e90: arXiv:2009.07803 Eq.(3.2) topology (c)              *)
    (* -------------------------------------------------------- *)
    "HB2L5P" -> <|
        "Description" -> "Non-Planar HexaBox 2L5P massless, s12=3,s23=2,s34=5,s45=7,s15=11",
        "Propagators" -> {
            l1^2,
            (l1 - p1)^2,
            (l1 - p1 - p2)^2,
            (l1 - p1 - p2 - p3)^2,
            l2^2,
            (l2 - p1 - p2 - p3 - p4)^2,
            (l1 - l2)^2,
            (l1 - l2 + p4)^2,
            (l2 - p1)^2,
            (l2 - p1 - p2)^2,
            (l2 - p1 - p2 - p3)^2
        },
        "LoopMomenta" -> {l1, l2},
        "ExternalMomenta" -> {p1, p2, p3, p4},
        "KinematicRules" -> {
            p1^2 -> 0, p2^2 -> 0, p3^2 -> 0, p4^2 -> 0,
            p1*p2 -> s12/2, p2*p3 -> s23/2, p3*p4 -> s34/2,
            p1*p3 -> (s45 - s12 - s23)/2,
            p2*p4 -> (s15 - s23 - s34)/2,
            p1*p4 -> (s23 - s45 - s15)/2
        },
        "TopSector" -> {1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0},
        "Numeric" -> {s12 -> 3, s23 -> 2, s34 -> 5, s45 -> 7, s15 -> 11, "d" -> 1/13},
        "Modulus" -> Prime[10000000],
        "Graph"->{{1,6},{2,1},{3,2},{7,3},{6,5},{5,7},{6,4},{4,7}}
    |>,

    (* -------------------------------------------------------- *)
    (* PB2L5P \[LongDash] Planar PentaBox 2L5P, massless                  *)
    (* \:6765\:6e90: arXiv:2009.07803 Eq.(3.2) topology (b)              *)
    (* -------------------------------------------------------- *)
    "PB2L5P" -> <|
        "Description" -> "Planar PentaBox 2L5P massless, s12=3,s23=2,s34=5,s45=7,s15=11",
        "Propagators" -> {
            l1^2,
            (l1 + p1)^2,
            (l1 + p1 + p2)^2,
            (l1 + p1 + p2 + p3)^2,
            l2^2,
            (l2 + p1 + p2 + p3)^2,
            (l2 + p1 + p2 + p3 + p4)^2,
            (l1 - l2)^2,
            (l1 + p1 + p2 + p3 + p4)^2,
            (l2 + p1)^2,
            (l2 + p1 + p2)^2
        },
        "LoopMomenta" -> {l1, l2},
        "ExternalMomenta" -> {p1, p2, p3, p4},
        "KinematicRules" -> {
            p1^2 -> 0, p2^2 -> 0, p3^2 -> 0, p4^2 -> 0,
            p1*p2 -> s12/2, p2*p3 -> s23/2, p3*p4 -> s34/2,
            p1*p3 -> (s45 - s12 - s23)/2,
            p2*p4 -> (s15 - s23 - s34)/2,
            p1*p4 -> (s23 - s45 - s15)/2
        },
        "TopSector" -> {1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0},
        "Numeric" -> {s12 -> 3, s23 -> 2, s34 -> 5, s45 -> 7, s15 -> 11, "d" -> 1/13},
        "Modulus" -> Prime[10000000],
        "Graph"->{{6,1},{1,2},{2,3},{3,7},{5,6},{7,4},{4,5},{7,6}}
    |>,

    (* -------------------------------------------------------- *)
    (* Penta1L \[LongDash] 1-Loop Pentagon 5pt, massless                  *)
    (* \:6765\:6e90: arXiv:2009.07803 Eq.(3.2) topology (a)              *)
    (* -------------------------------------------------------- *)
    "Penta1L" -> <|
        "Description" -> "1L Pentagon 5pt massless, s12=1,s23=2,s34=3,s45=5,s15=8",
        "Propagators" -> {
            l1^2,
            (l1 + p1)^2,
            (l1 + p1 + p2)^2,
            (l1 + p1 + p2 + p3)^2,
            (l1 + p1 + p2 + p3 + p4)^2
        },
        "LoopMomenta" -> {l1},
        "ExternalMomenta" -> {p1, p2, p3, p4},
        "KinematicRules" -> {
            p1^2 -> 0, p2^2 -> 0, p3^2 -> 0, p4^2 -> 0,
            p1*p2 -> s12/2, p2*p3 -> s23/2, p3*p4 -> s34/2,
            p1*p3 -> (s45 - s12 - s23)/2,
            p2*p4 -> (s15 - s23 - s34)/2,
            p1*p4 -> (s23 - s45 - s15)/2
        },
        "TopSector" -> {1, 1, 1, 1, 1},
        "Numeric" -> {s12 -> 1, s23 -> 2, s34 -> 3, s45 -> 5, s15 -> 8, "d" -> 1/3},
        "Modulus" -> Prime[10000000],
        "Graph"->{{5,1},{1,2},{2,3},{3,4},{4,5}}
    |>,

    (* -------------------------------------------------------- *)
    (* Penta1L-5m \[LongDash] 1L Pentagon 5pt, all-massive (msq=11)        *)
    (* -------------------------------------------------------- *)
    "Penta1L-5m" -> <|
        "Description" -> "1L Pentagon 5pt all-massive, s12=1,s23=2,s34=3,s45=5,s15=8,msq=11",
        "Propagators" -> {
            l1^2 - msq,
            (l1 + p1)^2 -  msq,
            (l1 + p1 + p2)^2 - msq,
            (l1 + p1 + p2 + p3)^2 - msq,
            (l1 + p1 + p2 + p3 + p4)^2 - msq
        },
        "LoopMomenta" -> {l1},
        "ExternalMomenta" -> {p1, p2, p3, p4},
        "KinematicRules" -> {
            p1^2 -> 0, p2^2 -> 0, p3^2 -> 0, p4^2 -> 0,
            p1*p2 -> s12/2, p2*p3 -> s23/2, p3*p4 -> s34/2,
            p1*p3 -> (s45 - s12 - s23)/2,
            p2*p4 -> (s15 - s23 - s34)/2,
            p1*p4 -> (s23 - s45 - s15)/2
        },
        "TopSector" -> {1, 1, 1, 1, 1},
        "Numeric" -> {s12 -> 1, s23 -> 2, s34 -> 3, s45 -> 5, s15 -> 8, msq -> 11, "d" -> 1/3},
        "Modulus" -> Prime[10000000],
        "Graph"->{{5,1},{1,2},{2,3},{3,4},{4,5}}
    |>,

    (* -------------------------------------------------------- *)
    (* LA3L4Pm \[LongDash] 3L Ladder A 4pt, 2-massive                      *)
    (* \:6765\:6e90: arXiv:2410.15431 Eq.(1) Family A                    *)
    (* -------------------------------------------------------- *)
    "LA3L4Pm" -> <|
        "Description" -> "3L Ladder A 4pt 2-massive, s=3,t=2,m=1",
        "Propagators" -> {
            l1^2,
            (l1 + p1)^2,
            (l1 + p1 + p2)^2,
            (l2 + p1 + p2)^2,
            (l3 + p1 + p2)^2,
            (l3 + p1 + p2 + p3)^2,
            l3^2,
            l2^2,
            (l1 - l2)^2,
            (l2 - l3)^2,
            (l1 + p1 + p2 + p3)^2,
            (l2 + p1 + p2 + p3)^2,
            (l2 + p1)^2,
            (l3 + p1)^2,
            (l1 - l3)^2
        },
        "LoopMomenta" -> {l1, l2, l3},
        "ExternalMomenta" -> {p1, p2, p3},
        "KinematicRules" -> {
            p1^2 -> 0, p2^2 -> 0, p3^2 -> msq,
            p1*p2 -> s/2,
            p2*p3 -> (t - msq)/2,
            p1*p3 -> (msq - s - t)/2
        },
        "TopSector" -> {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0},
        "Numeric" -> {s -> 3, t -> 2, msq -> 1, "d" -> 1/3},
        "Modulus" -> Prime[10000000],
        "Graph"->{{5,1},{1,2},{2,6},{6,7},{7,3},{3,4},{4,8},{8,5},{6,5},{7,8}}
    |>,

    (* -------------------------------------------------------- *)
    (* LB3L4Pm \[LongDash] 3L Ladder B 4pt, 2-massive                      *)
    (* \:6765\:6e90: arXiv:2410.15431 Eq.(2) Family B                    *)
    (* -------------------------------------------------------- *)
    "LB3L4Pm" -> <|
        "Description" -> "3L Ladder B 4pt 2-massive, s=3,t=2,m=1",
        "Propagators" -> {
            l1^2,
            (l1 + p1)^2,
            (l2 + p1)^2,
            (l3 + p1)^2,
            (l3 + p1 + p2)^2,
            (l3 + p1 + p2 + p3)^2,
            (l2 + p1 + p2 + p3)^2,
            (l1 + p1 + p2 + p3)^2,
            (l1 - l2)^2,
            (l2 - l3)^2,
            (l1 + p1 + p2)^2,
            (l2 + p1 + p2)^2,
            l3^2,
            l2^2,
            (l1 - l3)^2
        },
        "LoopMomenta" -> {l1, l2, l3},
        "ExternalMomenta" -> {p1, p2, p3},
        "KinematicRules" -> {
            p1^2 -> 0, p2^2 -> 0, p3^2 -> msq,
            p1*p2 -> s/2,
            p2*p3 -> (t - msq)/2,
            p1*p3 -> (msq - s - t)/2
        },
        "TopSector" -> {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0},
        "Numeric" -> {s -> 3, t -> 2, msq -> 1, "d" -> 1/3},
        "Modulus" -> Prime[10000000],
        "Graph"->{{4,1},{1,5},{5,6},{6,2},{2,3},{3,7},{7,8},{8,4},{5,8},{6,7}}
    |>

|>;

(* ------------------------------------------------------------ *)
(* JSON \:540c\:6b65                                                      *)
(* ------------------------------------------------------------ *)

(* \:5c06\:5355\:4e2a\:79ef\:5206\:65cf\:5bfc\:51fa\:5230 families/<Family>.json *)
ExportFamilyToJSON[family_String] := Module[
    {cfg, jsonData, jsonPath, familiesDir, kinRules, numericRules, propagators,
     loopMomenta, externalMomenta, topSector},
    cfg = $FamilyDatabase[family];
    familiesDir = If[$InputFileName =!= "",
        DirectoryName[$InputFileName],
        FileNameJoin[{Directory[], "families"}]
    ];

    (* \:4f20\:64ad\:5b50 -> \:5b57\:7b26\:4e32\:6570\:7ec4 *)
    propagators = ToString[#, InputForm] & /@ cfg["Propagators"];

    (* \:5708\:52a8\:91cf -> \:5b57\:7b26\:4e32\:6570\:7ec4 *)
    loopMomenta = ToString[#] & /@ cfg["LoopMomenta"];

    (* \:5916\:52a8\:91cf -> \:5b57\:7b26\:4e32\:6570\:7ec4 *)
    externalMomenta = ToString[#] & /@ cfg["ExternalMomenta"];

    (* \:8fd0\:52a8\:5b66\:89c4\:5219 -> <|"expr" -> "expr"|> *)
    kinRules = Association[
        (ToString[#[[1]], InputForm] -> ToString[#[[2]], InputForm]) & /@
         Normal[cfg["KinematicRules"]]
    ];

    (* \:6570\:503c\:89c4\:5219\:ff1a\:6574\:6570\:4fdd\:7559\:539f\:503c\:ff0c\:6709\:7406\:6570/\:5176\:4ed6\:8f6c\:5b57\:7b26\:4e32 *)
    numericRules = Association[
        Function[r,
            With[{k = If[StringQ[r[[1]]], r[[1]], ToString[r[[1]]]]},
                k -> Switch[r[[2]],
                    _Integer | _Real, r[[2]],
                    _, ToString[r[[2]], InputForm]
                ]
            ]
        ] /@ Normal[cfg["Numeric"]]
    ];

    (* \:9876\:6247\:533a *)
    topSector = cfg["TopSector"];

    jsonData = <|
        "name" -> family,
        "description" -> cfg["Description"],
        "propagators" -> propagators,
        "loopMomenta" -> loopMomenta,
        "externalMomenta" -> externalMomenta,
        "kinematicRules" -> kinRules,
        "topSector" -> topSector,
        "numeric" -> numericRules,
        "modulus" -> cfg["Modulus"],
        "graph" -> cfg["Graph"]
    |>;

    jsonPath = FileNameJoin[{familiesDir, family <> ".json"}];
    Export[jsonPath, jsonData, "JSON", "Compact" -> False];
    Print["  \[Checkmark] ", family, " \[RightArrow] ", family <> ".json"];
];

(* \:5bfc\:51fa\:6240\:6709\:79ef\:5206\:65cf\:5230 families/*.json *)
ExportAllFamiliesToJSON[] := Module[{fams},
    fams = ListFamilies[];
    Print["=== Exporting ", Length[fams], " families to JSON ==="];
    Scan[ExportFamilyToJSON, fams];
    Print["=== Done ==="];
];

(* ------------------------------------------------------------ *)
(* \:8f85\:52a9\:51fd\:6570                                                       *)
(* ------------------------------------------------------------ *)

(* \:5217\:51fa\:6240\:6709\:5df2\:6ce8\:518c\:7684\:79ef\:5206\:65cf *)
ListFamilies[] := Keys[$FamilyDatabase];

(* \:83b7\:53d6\:79ef\:5206\:65cf\:7684\:4f20\:64ad\:5b50\:4e2a\:6570 *)
GetNProp[family_String] := Module[{cfg},
    cfg = $FamilyDatabase[family];
    Length[cfg["Propagators"]]
];

(* \:83b7\:53d6\:79ef\:5206\:65cf\:7684\:5916\:817f\:6570 *)
GetNE[family_String] := Module[{cfg},
    cfg = $FamilyDatabase[family];
    Length[cfg["ExternalMomenta"]]
];

(* \:83b7\:53d6\:79ef\:5206\:65cf\:7684\:5708\:6570 *)
GetNLoop[family_String] := Module[{cfg},
    cfg = $FamilyDatabase[family];
    Length[cfg["LoopMomenta"]]
];

(* \:4e3a LIEDefineFamily \:5c55\:5f00\:914d\:7f6e *)
MakeFamilyConfig[family_String] := Module[{cfg, props},
    cfg = $FamilyDatabase[family];
    props = cfg["Propagators"] /. cfg["Numeric"];
    <|
        cfg,
        "Propagators" -> props,
        "KinematicRules" -> (cfg["KinematicRules"] /. cfg["Numeric"])
    |>
];

(* \:6253\:5370\:79ef\:5206\:65cf\:4fe1\:606f\:6458\:8981 *)
PrintFamilyInfo[family_String] := Module[{cfg},
    cfg = $FamilyDatabase[family];
    Print["Family: ", family];
    Print["  Description: ", cfg["Description"]];
    Print["  Loops: ", Length[cfg["LoopMomenta"]], " (", cfg["LoopMomenta"], ")"];
    Print["  External: ", Length[cfg["ExternalMomenta"]], " (", cfg["ExternalMomenta"], ")"];
    Print["  Propagators: ", Length[cfg["Propagators"]]];
    Print["  TopSector: ", cfg["TopSector"]];
    Print["  Numeric: ", cfg["Numeric"]];
    Print["  Modulus: ", cfg["Modulus"]];
];

(* \:6253\:5370\:6240\:6709\:5df2\:6ce8\:518c\:79ef\:5206\:65cf *)
PrintAllFamilies[] := Module[{fams},
    fams = ListFamilies[];
    Print["=== Registered IBP Families (", Length[fams], ") ==="];
    Do[Print["  ", f, " \[LongDash] ", $FamilyDatabase[f]["Description"]], {f, fams}];
];

(* ------------------------------------------------------------ *)
(* \:521d\:59cb\:5316\:8f93\:51fa                                                      *)
(* ------------------------------------------------------------ *)
PrintAllFamilies[];
