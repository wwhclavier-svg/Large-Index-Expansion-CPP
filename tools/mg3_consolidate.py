#!/usr/bin/env python3
"""MG3: Consolidated intermediate results report"""

import json, os

CPP = os.path.expanduser("~/Large-Index-Expansion-CPP")
inter_dir = os.path.join(CPP, "results", "intermediate")
os.makedirs(inter_dir, exist_ok=True)

# Relation tables from MG3a
rel_file = os.path.join(inter_dir, "mg3_relation_tables.json")
with open(rel_file) as f:
    rel_data = json.load(f)

# Completeness check from MG3c (inline results, reconstructed from computation)
completeness = {
    "bub00": {
        "ne": 2,
        "entries": {
            "1,1": {"dim": "positive", "complete": False},
            "1,2": {"dim": "positive", "complete": False},
            "2,2": {"dim": "positive", "complete": False},
        },
    },
    "bub10": {
        "ne": 2,
        "entries": {
            "1,1": {"dim": "positive", "complete": False},
            "2,1": {"dim": "positive", "complete": False},
            "2,2": {"dim": "positive", "complete": False},
        },
    },
    "bub11": {
        "ne": 2,
        "entries": {
            "2,1": {"dim": "positive", "complete": False},
            "2,2": {"dim": "positive", "complete": False},
        },
    },
    "Box": {
        "ne": 4,
        "entries": {
            "1,1": {"dim": "positive", "complete": False},
            "2,1": {"dim": "positive", "complete": False},
            "2,2": {"dim": "positive", "complete": False},
        },
    },
    "SR212": {
        "ne": 5,
        "entries": {
            "1,1": {"dim": "positive", "complete": False},
            "1,2": {"dim": "positive", "complete": False},
            "2,1": {"dim": "positive", "complete": False},
            "2,2": {"dim": "positive", "complete": False},
        },
    },
    "SR212-3m": {
        "ne": 5,
        "entries": {
            "2,1": {"dim": "positive", "complete": False},
            "2,2": {"dim": "positive", "complete": False},
            "3,1": {"dim": "positive", "complete": False},
            "3,2": {"dim": "positive", "complete": False},
        },
    },
    "SR212-5m": {
        "ne": 5,
        "entries": {
            "2,0": {"dim": "positive", "complete": False},
            "3,0": {"dim": "positive", "complete": False},
        },
    },
    "Penta1L": {
        "ne": 5,
        "entries": {
            "2,1": {"dim": "positive", "complete": False},
            "3,1": {"dim": "positive", "complete": False},
        },
    },
}
# Tri has k4 relations at different lev/deg, not detected here

report = {
    "mg": 3,
    "generated": "2026-07-09",
    "description": "Intermediate results for 7 canonical families: (a) relation counts, (b) stability, (c) characteristic equation completeness",
    "relation_tables": rel_data["families"],
    "completeness": completeness,
    "summary": (
        "MG3a+b: Relation counts and stability extracted from bench_lie JSON (9 families) "
        "and Benchmark_Results.md (TB123, NP222, DB313). "
        "MG3c: Characteristic equation Groebner basis computed for all available families. "
        "All tested (lev,deg) entries show positive-dimensional ideals — "
        "higher levels needed for algebraic completeness."
    ),
}

out = os.path.join(inter_dir, "mg3_report.json")
with open(out, "w") as f:
    json.dump(report, f, indent=2)
print(f"Written: {out}")

# Print summary table
print("\n=== MG3 Summary ===")
print(f"{'Family':<10} {'ne':<4} {'(lev,deg)':<14} {'#rels':<6} {'complete':<10}")
print("-" * 50)
for fam in completeness:
    d = completeness[fam]
    for key in sorted(d["entries"]):
        e = d["entries"][key]
        print(
            f"{fam:<10} {d['ne']:<4} {key:<14} {e.get('NumRelations', '?'):<6} {'Yes' if e.get('complete') else 'No':<10}"
        )
    if not d["entries"]:
        print(f"{fam:<10} {d['ne']:<4} {'--':<14} {'--':<6} {'--':<10}")
