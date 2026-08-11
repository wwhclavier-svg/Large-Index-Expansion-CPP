#!/usr/bin/env python3
"""MG3a: Consolidate relation counts + stable orders from bench_lie JSON + Benchmark_Results.md"""

import json, re, os, glob

CPP = os.path.expanduser("~/Large-Index-Expansion-CPP")
lie_dir = os.path.join(CPP, "results", "lie")
out_dir = os.path.join(CPP, "results", "intermediate")
os.makedirs(out_dir, exist_ok=True)

# ============================================================
# 1. Extract from bench_lie JSONs
# ============================================================
json_files = glob.glob(os.path.join(lie_dir, "*.json"))
json_data = {}
for f in sorted(json_files):
    with open(f) as fh:
        j = json.load(fh)
    fam = j["family"]
    counts = j.get("relation_counts", {})
    levels = j.get("phases", {}).get("relation_solve", {}).get("levels", [])
    # Build stability table from levels
    stability = {}
    for lev in levels:
        l = lev["lev"]
        d = lev["deg"]
        stability.setdefault(l, {})[d] = {
            "solutions": lev.get("solution_dimension", 0),
            "stable_order": lev.get("stable_order", 0),
            "converged": lev.get("converged", False),
            "active_vars": lev.get("active_variables", 0),
            "sampling_points": lev.get("sampling_points_used", 0),
        }
    json_data[fam] = {
        "ne": j.get("config", {}).get("ne", 0),
        "order": j.get("config", {}).get("order", 0),
        "prime": j.get("config", {}).get("prime", 0),
        "relation_counts": counts.get("counts", {}),
        "stability": stability,
    }
    print(f"[JSON] {fam}: ne={json_data[fam]['ne']}")

# ============================================================
# 2. Parse Benchmark_Results.md for TB123, NP222, DB313
# ============================================================
md_file = os.path.join(CPP, "docs", "reports", "Benchmark_Results.md")
with open(md_file) as fh:
    md_text = fh.read()

# Parse: find sections with relation counts and stable orders
# The format is like:
# | (lev,deg) | N | stable_order |
# We extract known values from the existing doc
md_families = {
    "TB123": {
        "ne": 7,
        "sections": {
            "topsector k=4": {
                "order": 4,
                "mode": "topsector",
                "counts": {"lev3": {"deg0": 21}},
                "stability": {"lev3": {"deg0": {"stable_order": 0, "converged": True}}},
            },
            "full k=2": {
                "order": 2,
                "mode": "full",
                "counts": {"lev3": {"deg0": 21, "deg2": 284}},
                "stability": {
                    "lev3": {
                        "deg0": {"stable_order": 0, "converged": True},
                        "deg2": {"stable_order": -2, "converged": False},
                    }
                },
            },
        },
    },
    "NP222": {
        "ne": 7,
        "sections": {
            "topsector k=4": {
                "order": 4,
                "mode": "topsector",
                "counts": {"lev3": {"deg0": 21, "deg1": 392, "deg2": 1862}},
                "stability": {
                    "lev3": {
                        "deg0": {"stable_order": 1, "converged": True},
                        "deg1": {"stable_order": 3, "converged": True},
                        "deg2": {"stable_order": -2, "converged": False},
                    }
                },
            }
        },
    },
    "DB313": {
        "ne": 9,
        "sections": {
            "topsector k=2": {
                "order": 2,
                "mode": "topsector",
                "counts": {"lev3": {"deg0": 36, "deg1": 1168}},
                "stability": {
                    "lev3": {
                        "deg0": {"stable_order": 0, "converged": True},
                        "deg1": {"stable_order": -2, "converged": False},
                    }
                },
            },
            "full k=2": {
                "order": 2,
                "mode": "full",
                "counts": {"lev3": {"deg0": 36, "deg1": 0}},
                "stability": {
                    "lev3": {
                        "deg0": {"stable_order": 0, "converged": True},
                        "deg1": {"stable_order": 1, "converged": True},
                    }
                },
            },
            "topsector k=3": {
                "order": 3,
                "mode": "topsector",
                "counts": {"lev3": {"deg0": 36, "deg1": 568}},
                "stability": {
                    "lev3": {
                        "deg0": {"stable_order": 0, "converged": True},
                        "deg1": {"stable_order": -2, "converged": False},
                    }
                },
            },
            "full k=3": {
                "order": 3,
                "mode": "full",
                "counts": {"lev3": {"deg0": 36}},
                "stability": {"lev3": {"deg0": {"stable_order": 1, "converged": True}}},
            },
        },
    },
}

# Also add Penta1L-5m from JSON
# (already in json_data from the bench_lie scan)

# ============================================================
# 3. Build unified output
# ============================================================
all_results = {}

# First add JSON-derived families
for fam, data in sorted(json_data.items()):
    all_results[fam] = {
        "ne": data["ne"],
        "order": data["order"],
        "prime": data["prime"],
        "source": "bench_lie",
        "relation_counts": data["relation_counts"],
        "stability": data["stability"],
    }

# Add markdown-derived families
for fam, info in md_families.items():
    sections = {}
    for name, sec in info["sections"].items():
        sections[name] = {
            "order": sec["order"],
            "mode": sec["mode"],
            "relation_counts": sec["counts"],
            "stability": sec["stability"],
        }
    all_results[fam] = {
        "ne": info["ne"],
        "source": "Benchmark_Results.md",
        "sections": sections,
    }

# Write consolidated JSON
out_file = os.path.join(out_dir, "mg3_relation_tables.json")
with open(out_file, "w") as fh:
    json.dump({"generated": "2026-07-09", "families": all_results}, fh, indent=2)
print(f"\nWrote: {out_file}")

# Print summary table
print("\n=== Relation Counts Summary ===")
print(f"{'Family':<12} {'ne':<4} {'Config':<20} {'Counts'}")
print("-" * 60)
for fam in sorted(all_results.keys()):
    d = all_results[fam]
    if d["source"] == "bench_lie":
        cnt = d["relation_counts"]
        entries = []
        for l in sorted(cnt.keys()):
            for dd in sorted(cnt[l].keys()):
                if cnt[l][dd] > 0:
                    entries.append(f"({l[-1]},{dd})={cnt[l][dd]}")
        print(
            f"{fam:<12} {d['ne']:<4} {'order=' + str(d['order']):<20} {', '.join(entries)}"
        )
    else:
        for sname, sec in d["sections"].items():
            cnt = sec["relation_counts"]
            entries = []
            for l in sorted(cnt.keys()):
                for dd in sorted(cnt[l].keys()):
                    if cnt[l][dd] > 0:
                        entries.append(f"({l[-1]},{dd})={cnt[l][dd]}")
            print(f"{fam:<12} {d['ne']:<4} {sname:<20} {', '.join(entries)}")
