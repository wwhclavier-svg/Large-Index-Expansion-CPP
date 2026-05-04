(* ============================================================ *)
(* FamilyDatabase.wl — IBP 积分族集中定义数据库                    *)
(* ============================================================ *)
(* 用法:                                                           *)
(*   Get["verify/FamilyDatabase/FamilyDatabase.wl"];               *)
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
    (* 1-Loop 族 (L=1)                                          *)
    (* ======================================================== *)

    (* -------------------------------------------------------- *)
    (* bub00 — 1-loop bubble, massless (msq=0)                  *)
    (* 传播子: {-k1^2, -(k1-p1)^2}                               *)
    (* -------------------------------------------------------- *)
    "bub00" -> <|
        "Description" -> "1-loop bubble, msq=0",
        "Propagators" -> {-k1^2, -(k1 - p1)^2},
        "LoopMomenta" -> {k1},
        "ExternalMomenta" -> {p1},
        "KinematicRules" -> {p1^2 -> s},
        "TopSector" -> {1, 1},
        "Numeric" -> {s -> 3, msq -> 0, "d" -> 1/3},
        "Modulus" -> Prime[10000000]  (* 179424673 *)
    |>,

    (* -------------------------------------------------------- *)
    (* bub10 — 1-loop bubble, 1-massive (msq=1)                 *)
    (* 传播子: {-k1^2+msq, -(k1-p1)^2}                           *)
    (* -------------------------------------------------------- *)
    "bub10" -> <|
        "Description" -> "1-loop bubble, msq=1",
        "Propagators" -> {-k1^2 + msq, -(k1 - p1)^2},
        "LoopMomenta" -> {k1},
        "ExternalMomenta" -> {p1},
        "KinematicRules" -> {p1^2 -> s},
        "TopSector" -> {1, 1},
        "Numeric" -> {s -> 3, msq -> 1, "d" -> 1/3},
        "Modulus" -> Prime[10000000]
    |>,

    (* -------------------------------------------------------- *)
    (* bub11 — 1-loop bubble, both-massive (msq=1)              *)
    (* 传播子: {-k1^2 + msq, -(k1-p1)^2 + msq}                   *)
    (* -------------------------------------------------------- *)
    "bub11" -> <|
        "Description" -> "1-loop bubble, msq=1",
        "Propagators" -> {-k1^2 + msq, -(k1 - p1)^2 + msq},
        "LoopMomenta" -> {k1},
        "ExternalMomenta" -> {p1},
        "KinematicRules" -> {p1^2 -> s},
        "TopSector" -> {1, 1},
        "Numeric" -> {s -> 3, msq -> 1, "d" -> 1/3},
        "Modulus" -> Prime[10000000]
    |>,

    (* -------------------------------------------------------- *)
    (* Tri — 1-loop triangle (3 propagators, 2 ext momenta)     *)
    (* Kinematics: p1^2=s1, p2^2=s2, (p1+p2)^2=s3               *)
    (* -------------------------------------------------------- *)
    "Tri" -> <|
        "Description" -> "1-loop triangle, 3 props, s1=3,s2=8,s3=0",
        "Propagators" -> {-k1^2 + msq, -(k1 - p1)^2 + msq, -(k1 + p2)^2 + msq},
        "LoopMomenta" -> {k1},
        "ExternalMomenta" -> {p1, p2},
        "KinematicRules" -> {
            p1^2 -> s1, p2^2 -> s2,
            (p1 + p2)^2 -> s3,
            p1*p2 -> (s3 - s1 - s2)/2
        },
        "TopSector" -> {1, 1, 1},
        "Numeric" -> {s1 -> 3, s2 -> 8, s3 -> 0, msq -> 0, "d" -> 1/3},
        "Modulus" -> Prime[10000000]
    |>,

    (* -------------------------------------------------------- *)
    (* Box — 1-loop box (4 propagators, 3 ext momenta)          *)
    (* Kinematics: p1^2=s, p2^2=t, p3^2=u, Mandelstam          *)
    (* -------------------------------------------------------- *)
    "Box" -> <|
        "Description" -> "1-loop box, 4 props, s=3,t=2,u=-5",
        "Propagators" -> {
            -k1^2 + msq,
            -(k1 - p1)^2 + msq,
            -(k1 - p1 - p2)^2 + msq,
            -(k1 - p1 - p2 - p3)^2 + msq
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
        "Modulus" -> Prime[10000000]
    |>,

    (* ======================================================== *)
    (* 2-Loop 族 (L=2)                                          *)
    (* ======================================================== *)

    (* -------------------------------------------------------- *)
    (* SR — Sunrise (2L1P, 5 propagators)                       *)
    (* 来源: LIETest_2L2P_SR212.nb "Massless" 节                *)
    (* -------------------------------------------------------- *)
    "SR" -> <|
        "Description" -> "Sunrise 2L1P, massless, s=1",
        "Propagators" -> {
            -l1^2,
            -(l1 + p)^2,
            -l2^2,
            -(l2 + p)^2,
            -(l1 + l2 + p)^2
        },
        "LoopMomenta" -> {l1, l2},
        "ExternalMomenta" -> {p},
        "KinematicRules" -> {p^2 -> s},
        "TopSector" -> {1, 1, 1, 1, 1},
        "Numeric" -> {s -> 1, "d" -> 1/13},
        "Modulus" -> Prime[10000000]
    |>,
    "SR212" -> <|
        "Description" -> "Sunrise 2L1P, massless, s=1",
        "Propagators" -> {
            -l1^2,
            -(l1 + p)^2,
            -l2^2,
            -(l2 + p)^2,
            -(l1 + l2 + p)^2
        },
        "LoopMomenta" -> {l1, l2},
        "ExternalMomenta" -> {p},
        "KinematicRules" -> {p^2 -> s},
        "TopSector" -> {1, 1, 1, 1, 1},
        "Numeric" -> {s -> 1, "d" -> 1/13},
        "Modulus" -> Prime[10000000]
    |>,

    (* -------------------------------------------------------- *)
    "SR3m" -> <|
        "Description" -> "Sunrise 2L1P, 3-massive (props 1,3,5), s=3/2",
        "Propagators" -> {
            -l1^2 - msq,
            -(l1 + p)^2,
            -l2^2 - msq,
            -(l2 + p)^2,
            -(l1 + l2 + p)^2 - msq
        },
        "LoopMomenta" -> {l1, l2},
        "ExternalMomenta" -> {p},
        "KinematicRules" -> {p^2 -> s},
        "TopSector" -> {1, 0, 1, 0, 1},
        "Numeric" -> {s -> 3/2, msq -> 1, "d" -> 1/7},
        "Modulus" -> Prime[10000000]
    |>,
    (* SR212-3m — Sunrise 3-Massive (2L1P, 5 prop)                  *)
    (* 来源: LIETest_2L2P_SR212_3Massive.nb                     *)
    (* -------------------------------------------------------- *)
    "SR212-3m" -> <|
        "Description" -> "Sunrise 2L1P, 3-massive (props 1,3,5), s=3/2",
        "Propagators" -> {
            -l1^2 - msq,
            -(l1 + p)^2,
            -l2^2 - msq,
            -(l2 + p)^2,
            -(l1 + l2 + p)^2 - msq
        },
        "LoopMomenta" -> {l1, l2},
        "ExternalMomenta" -> {p},
        "KinematicRules" -> {p^2 -> s},
        "TopSector" -> {1, 0, 1, 0, 1},
        "Numeric" -> {s -> 3/2, msq -> 1, "d" -> 1/7},
        "Modulus" -> Prime[10000000]
    |>,

    (* -------------------------------------------------------- *)
    (* SR212-5m — Sunrise All-Massive (2L1P, 5 prop)                *)
    (* 来源: LIETest_2L2P_SR212.nb "All Massive" 节             *)
    (* -------------------------------------------------------- *)
    "SR212-5m" -> <|
        "Description" -> "Sunrise 2L1P, all-massive, s=0",
        "Propagators" -> {
            -l1^2 - msq,
            -(l1 + p)^2 - msq,
            -l2^2 - msq,
            -(l2 + p)^2 - msq,
            -(l1 + l2 + p)^2 - msq
        },
        "LoopMomenta" -> {l1, l2},
        "ExternalMomenta" -> {p},
        "KinematicRules" -> {p^2 -> s},
        "TopSector" -> {1, 1, 1, 1, 1},
        "Numeric" -> {s -> 0, msq -> 1, "d" -> 1/13},
        "Modulus" -> Prime[10000000]
    |>,

    (* -------------------------------------------------------- *)
    (* NP222 — Non-Planar 2-loop 2-point (2L2P, 7 propagators)  *)
    (* 来源: LIETest_2L2P_SR212.nb "Non Planar 2-loop" 节      *)
    (* -------------------------------------------------------- *)
    "NP222" -> <|
        "Description" -> "Non-Planar 2L2P, s=1,s2=0,m=0",
        "Propagators" -> {
            -(l2 + p1)^2 + msq,
            -(l1 - l2 - p1 - p2)^2,
            -l2^2 + msq,
            -(l1 - p2)^2 + msq,
            -(l1 - l2)^2,
            -l1^2 + msq,
            -(l1 + p1)^2
        },
        "LoopMomenta" -> {l1, l2},
        "ExternalMomenta" -> {p1, p2},
        "KinematicRules" -> {
            p1^2 -> 0, p2^2 -> s2,
            p1*p2 -> (s - s2)/2
        },
        "TopSector" -> {1, 1, 1, 1, 1, 1, 1},
        "Numeric" -> {s -> 1, s2 -> 0, msq -> 0, "d" -> 1/13},
        "Modulus" -> Prime[10000000]
    |>,

    (* -------------------------------------------------------- *)
    (* TB123 — Tri-Box / Tennis Ball 2L3P (7 prop, 6主+1辅)     *)
    (* 来源: LIETest_2L3P_TB123&NP222_Eigenvalue.nb Cell 3,5   *)
    (* -------------------------------------------------------- *)
    "TB123" -> <|
        "Description" -> "Tri-Box 2L3P, massless, s=3,s1=2,s2=4",
        "Propagators" -> {
            -l1^2,                         (* 1 *)
            -(l1 - p1)^2,                  (* 2 *)
            -(l1 - p1 - p2)^2,             (* 3 *)
            -(l2 + p1 + p2)^2,             (* 4 *)
            -l2^2,                         (* 5 *)
            -(l1 + l2)^2,                  (* 6 *)
            -(l2 + p1)^2                   (* 7: 辅助 ISP *)
        },
        "LoopMomenta" -> {l1, l2},
        "ExternalMomenta" -> {p1, p2},
        "KinematicRules" -> {
            p1^2 -> s1, p2^2 -> s2,
            p1*p2 -> (s - s1 - s2)/2
        },
        "TopSector" -> {1, 1, 1, 1, 1, 1, 0},
        "Numeric" -> {s -> 3, s1 -> 2, s2 -> 4, msq -> 0, "d" -> 1/13},
        "Modulus" -> Prime[10000000]
    |>,

    (* -------------------------------------------------------- *)
    (* TB123m — Tri-Box / Tennis Ball 2L3P massive              *)
    (* 来源: LIETest_2L3P_TB123&NP222_Eigenvalue.nb Cell 9,11  *)
    (* -------------------------------------------------------- *)
    "TB123m" -> <|
        "Description" -> "Tri-Box 2L3P, massive, s=3,s1=2,s2=5,m=1",
        "Propagators" -> {
            -l1^2 - msq,
            -(l1 - p1)^2 - msq,
            -(l1 - p1 - p2)^2 - msq,
            -(l2 + p1 + p2)^2 - msq,
            -l2^2 - msq,
            -(l1 + l2)^2 - msq,
            -(l2 + p1)^2
        },
        "LoopMomenta" -> {l1, l2},
        "ExternalMomenta" -> {p1, p2},
        "KinematicRules" -> {
            p1^2 -> s1, p2^2 -> s2,
            p1*p2 -> (s - s1 - s2)/2
        },
        "TopSector" -> {1, 1, 1, 1, 1, 1, 0},
        "Numeric" -> {s -> 3, s1 -> 2, s2 -> 5, msq -> 1, "d" -> 1/13},
        "Modulus" -> Prime[10000000]
    |>,

    (* -------------------------------------------------------- *)
    (* DB313 — Double Box (2L4P, 9 prop, 7主+2辅)              *)
    (* 来源: LIETest_2L4P_DB313_New.nb, dbox_data.wl           *)
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
            (l2 + p2)^2,            (* 8: 辅助 *)
            (l1 + p4)^2             (* 9: 辅助 *)
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
        "Modulus" -> Prime[10000000]
    |>,

    (* -------------------------------------------------------- *)
    (* NP322 — Non-Planar 2L4P (9 propagators)                  *)
    (* 来源: LIETest_2L4P_NP322.nb Cell 1-3                     *)
    (* 与 DB313 (Double Box) 不同: NP322 是非平面对应物          *)
    (* -------------------------------------------------------- *)
    "NP322" -> <|
        "Description" -> "Non-Planar 2L4P, massless, s=1,t=5",
        "Propagators" -> {
            -(k2 + p1)^2 + msq,             (* 1 *)
            -(k1 - k2 - p1 - p2)^2,         (* 2 *)
            -k2^2 + msq,                    (* 3 *)
            -(k1 - p2)^2 + msq,             (* 4 *)
            -(k1 - k2)^2,                   (* 5 *)
            -k1^2 + msq,                    (* 6 *)
            -(k1 - k2 + p4)^2,              (* 7 *)
            -(k1 + p1)^2,                   (* 8 *)
            -(k1 + p4)^2                    (* 9 *)
        },
        "LoopMomenta" -> {k1, k2},
        "ExternalMomenta" -> {p1, p2, p4},
        "KinematicRules" -> {
            p1^2 -> m1, p2^2 -> m2, p4^2 -> 0,
            p1*p2 -> (s - m1 - m2)/2,
            p1*p4 -> (t - m1)/2,
            p2*p4 -> -(s + t - m1)/2
        },
        "TopSector" -> {1, 1, 1, 1, 1, 1, 1, 1, 1},
        "Numeric" -> {s -> 1, t -> 5, m1 -> 0, m2 -> 0, msq -> 0, "d" -> 1/7},
        "Modulus" -> Prime[10000000]
    |>,

    (* -------------------------------------------------------- *)
    (* NP322m — Non-Planar 2L4P massive variant                  *)
    (* 来源: LIETest_2L4P_NP322.nb Cell 1 (massive numeric)     *)
    (* -------------------------------------------------------- *)
    "NP322m" -> <|
        "Description" -> "Non-Planar 2L4P, massive, s=3,t=2,m1=1,m2=4,m=1",
        "Propagators" -> {
            -(k2 + p1)^2 + msq,
            -(k1 - k2 - p1 - p2)^2,
            -k2^2 + msq,
            -(k1 - p2)^2 + msq,
            -(k1 - k2)^2,
            -k1^2 + msq,
            -(k1 - k2 + p4)^2,
            -(k1 + p1)^2,
            -(k1 + p4)^2
        },
        "LoopMomenta" -> {k1, k2},
        "ExternalMomenta" -> {p1, p2, p4},
        "KinematicRules" -> {
            p1^2 -> m1, p2^2 -> m2, p4^2 -> 0,
            p1*p2 -> (s - m1 - m2)/2,
            p1*p4 -> (t - m1)/2,
            p2*p4 -> -(s + t - m1)/2
        },
        "TopSector" -> {1, 1, 1, 1, 1, 1, 1, 1, 1},
        "Numeric" -> {s -> 3, t -> 2, m1 -> 1, m2 -> 4, msq -> 1, "d" -> 1/3},
        "Modulus" -> Prime[10000000]
    |>,

    (* -------------------------------------------------------- *)
    (* DP323 — Double Pentagon (2L5P, 11 prop, 8主+3辅)        *)
    (* 来源: dpen_data.wl, LIETest_2L5P_DP323.nb                *)
    (* -------------------------------------------------------- *)
    "DP323" -> <|
        "Description" -> "Double Pentagon 2L5P, s12=4,s13=7,s14=8,s23=1,s24=2,m=3",
        "Propagators" -> {
            -k1^2 + msq,                            (* 1  *)
            -(k1 - p1)^2 + msq,                     (* 2  *)
            -(k1 - p1 - p2)^2 + msq,                (* 3  *)
            -k2^2 + msq,                            (* 4  *)
            -(k2 - p1 - p2 - p3)^2 + msq,           (* 5  *)
            -(k2 - p1 - p2 - p3 - p4)^2 + msq,      (* 6  *)
            -(k1 - k2)^2 + msq,                     (* 7  *)
            -(k1 - k2 + p3)^2 + msq,                (* 8  *)
            -(k1 - p1 - p2 - p3 - p4)^2,            (* 9: 辅助 ISP *)
            -(k2 - p1)^2,                           (* 10: 辅助 ISP *)
            -(k2 - p1 - p2)^2                       (* 11: 辅助 ISP *)
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
        "Modulus" -> Prime[10000000]
    |>,

    (* ======================================================== *)
    (* 3-Loop 族 (L=3)                                          *)
    (* 来源: LIETest_2L4P_TP412&NP322_Eigenvalue.nb             *)
    (* ======================================================== *)

    (* -------------------------------------------------------- *)
    (* BN3L — 3-Loop Banana vacuum (3L0P, 6 propagators)        *)
    (* 来源: TP412&NP322_Eigenvalue.nb Cell 50                   *)
    (* -------------------------------------------------------- *)
    "BN3L" -> <|
        "Description" -> "3-Loop Banana, vacuum, msq=1",
        "Propagators" -> {
            -l1^2 + msq,
            -l2^2 + msq,
            -l3^2 + msq,
            -(l1 - l2)^2 + msq,
            -(l2 - l3)^2 + msq,
            -(l3 - l1)^2 + msq
        },
        "LoopMomenta" -> {l1, l2, l3},
        "ExternalMomenta" -> {},
        "KinematicRules" -> {},
        "TopSector" -> {1, 1, 1, 1, 1, 1},
        "Numeric" -> {msq -> 1, "d" -> 1/13},
        "Modulus" -> Prime[10000000]
    |>,

    (* -------------------------------------------------------- *)
    (* BN3L1P — 3-Loop Banana 1-Point (9 prop, 8主+1辅)        *)
    (* 来源: TP412&NP322_Eigenvalue.nb Cell 56                   *)
    (* -------------------------------------------------------- *)
    "BN3L1P" -> <|
        "Description" -> "3-Loop Banana 1-Point, massless, s=1",
        "Propagators" -> {
            -l1^2,                 (* 1 *)
            -(l1 + p)^2,           (* 2 *)
            -(l2 + p)^2,           (* 3 *)
            -l2^2,                 (* 4 *)
            -l3^2,                 (* 5 *)
            -(l1 - l2)^2,          (* 6 *)
            -(l2 - l3)^2,          (* 7 *)
            -(l3 - l1)^2,          (* 8 *)
            -(l3 - p)^2            (* 9: 辅助 ISP *)
        },
        "LoopMomenta" -> {l1, l2, l3},
        "ExternalMomenta" -> {p},
        "KinematicRules" -> {p^2 -> s},
        "TopSector" -> {1, 1, 1, 1, 1, 1, 1, 1, 0},
        "Numeric" -> {s -> 1, "d" -> 1/13},
        "Modulus" -> Prime[10000000]
    |>,

    (* -------------------------------------------------------- *)
    (* BN3L1Pm — 3-Loop Banana 1-Point massive (9 prop, 8主+1辅) *)
    (* 来源: TP412&NP322_Eigenvalue.nb Cell 62                   *)
    (* -------------------------------------------------------- *)
    "BN3L1Pm" -> <|
        "Description" -> "3-Loop Banana 1-Point, massive, s=1,msq=1",
        "Propagators" -> {
            -l1^2 + msq,
            -(l1 + p)^2 + msq,
            -(l2 + p)^2 + msq,
            -l2^2 + msq,
            -l3^2 + msq,
            -(l1 - l2)^2 + msq,
            -(l2 - l3)^2 + msq,
            -(l3 - l1)^2 + msq,
            -(l3 - p)^2
        },
        "LoopMomenta" -> {l1, l2, l3},
        "ExternalMomenta" -> {p},
        "KinematicRules" -> {p^2 -> s},
        "TopSector" -> {1, 1, 1, 1, 1, 1, 1, 1, 0},
        "Numeric" -> {s -> 1, msq -> 1, "d" -> 1/13},
        "Modulus" -> Prime[10000000]
    |>
|>;

(* ------------------------------------------------------------ *)
(* 辅助函数                                                       *)
(* ------------------------------------------------------------ *)

(* 列出所有已注册的积分族 *)
ListFamilies[] := Keys[$FamilyDatabase];

(* 获取积分族的传播子个数 *)
GetNProp[family_String] := Module[{cfg},
    cfg = $FamilyDatabase[family];
    Length[cfg["Propagators"]]
];

(* 获取积分族的外腿数 *)
GetNE[family_String] := Module[{cfg},
    cfg = $FamilyDatabase[family];
    Length[cfg["ExternalMomenta"]]
];

(* 获取积分族的圈数 *)
GetNLoop[family_String] := Module[{cfg},
    cfg = $FamilyDatabase[family];
    Length[cfg["LoopMomenta"]]
];

(* 为 LIEDefineFamily 展开配置 *)
MakeFamilyConfig[family_String] := Module[{cfg, props},
    cfg = $FamilyDatabase[family];
    props = cfg["Propagators"] /. cfg["Numeric"];
    <|
        cfg,
        "Propagators" -> props,
        "KinematicRules" -> (cfg["KinematicRules"] /. cfg["Numeric"])
    |>
];

(* 打印积分族信息摘要 *)
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

(* 打印所有已注册积分族 *)
PrintAllFamilies[] := Module[{fams},
    fams = ListFamilies[];
    Print["=== Registered IBP Families (", Length[fams], ") ==="];
    Do[Print["  ", f, " — ", $FamilyDatabase[f]["Description"]], {f, fams}];
];

(* ------------------------------------------------------------ *)
(* 初始化输出                                                      *)
(* ------------------------------------------------------------ *)
PrintAllFamilies[];
