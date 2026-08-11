#!/usr/bin/env python3
"""
Kira input generator for LIE integral families.

Converts families/*.json -> Kira project directory with:
  config/integralfamilies.yaml
  config/kinematics.yaml
  jobs.yaml
  integrals

Usage:
  ./tools/kira_input_gen.py <family_name> [--output-dir kira_tests/]
  ./tools/kira_input_gen.py --all-1l2p     # all 1L + 2L2P families
"""

import json, os, sys, re

FAMILIES_DIR = os.path.join(os.path.dirname(__file__), "..", "families")
DEFAULT_OUTPUT = os.path.join(os.path.dirname(__file__), "..", "kira_tests")

ONELOOP = ["bub00", "bub10", "bub11", "Tri", "Box", "Penta1L", "Penta1L-5m"]
TWOLOOP_2POINT = ["SR212", "SR212-3m", "SR212-5m"]
ALL_1L2P = ONELOOP + TWOLOOP_2POINT


# ---------------------------------------------------------------------------
# Propagator conversion
# ---------------------------------------------------------------------------
def parse_propagator(prop_str, mass_symbol=None):
    s = prop_str.replace(" ", "")
    mass = 0
    # 有质量传播子 (msq≠0): Kira 约定 D = 表达式 - m^2。族定义给出 "P ± msq" (P=动量部分),
    # 为维持与无质量族一致的 (-1)^Σ 约定, 需整体取负: -(P ± msq) = -P ∓ msq,
    # 即 kira_expr = -P + (∓msq)。
    # 实测: ["-k1^2+msq", 0] 正确 (bub10 IBP 核验通过); 原实现 ["-k1^2","msq"]
    # 使 Kira 算成 D = -k1^2 - msq (质量符号错, 得到错误的族)。
    if "msq" in s and mass_symbol:
        if "-msq" in s:
            mass_sign = -1
            M = s.replace("-msq", "")
        else:
            mass_sign = 1
            M = s.replace("+msq", "").replace("msq", "")
        M = M.strip("+-")
        if M.startswith("-"):
            negM = M[1:]
        else:
            negM = "-" + M
        if mass_sign == -1:
            kira_expr = negM + "+" + mass_symbol
        else:
            kira_expr = negM + "-" + mass_symbol
        return [kira_expr, 0]
    # msq=0 (无质量): 剥离 msq 后整体取负, 维持 (-1)^Σ 约定
    if "msq" in s:
        s = s.replace("-msq", "").replace("+msq", "").replace("msq", "")
        s = s.strip("+-")
    if s.startswith("-"):
        kira_expr = s[1:]
    else:
        kira_expr = "-" + s
    return [kira_expr, mass]


# ---------------------------------------------------------------------------
# Kinematics helpers
# ---------------------------------------------------------------------------
def extract_invariants(kinematic_rules):
    invars = set()
    for val in kinematic_rules.values():
        for token in re.split(r"[+\-*/()\s]", val):
            token = token.strip()
            if token and token[0].isalpha() and token != "0":
                invars.add(token)
    return sorted(invars)


def scalar_to_kira_rule(key, val):
    lhs = key.replace(" ", "")
    if lhs.endswith("^2"):
        inner = lhs[:-2]
        if inner.startswith("(") and inner.endswith(")"):
            return None
        return [[inner, inner], val]
    elif "*" in lhs:
        parts = lhs.split("*")
        return [[parts[0], parts[1]], val]
    return None


def sector_to_int(sector_vec):
    # Kira 约定: bit i (2^i) 对应第 i+1 个传播子, 即列表首位是最低位
    # (此前未反转, 回文 sector 如 SR212 的 31 恰好正确, 非回文则错位)
    return int("".join(str(b) for b in reversed(sector_vec)), 2)


# ---------------------------------------------------------------------------
# YAML emitter (produces Kira-compatible format, properly indented)
# ---------------------------------------------------------------------------
def yaml_str(obj):
    """Convert Python object to Kira-compatible YAML string."""
    lines = []

    def emit(obj, indent=0, prefix=""):
        pad = "  " * indent
        if isinstance(obj, dict):
            for k, v in obj.items():
                if isinstance(v, list):
                    if not v:
                        lines.append(f"{pad}{k}: []")
                    elif isinstance(v[0], list):
                        lines.append(f"{pad}{k}:")
                        for item in v:
                            inner = ", ".join(q(x) for x in item)
                            lines.append(f"{pad}  - [{inner}]")
                    elif isinstance(v[0], dict):
                        lines.append(f"{pad}{k}:")
                        for item in v:
                            lines.append(f"{pad}  -")
                            emit(item, indent + 2)
                    else:
                        inner = ", ".join(q(x) for x in v)
                        lines.append(f"{pad}{k}: [{inner}]")
                elif isinstance(v, dict):
                    lines.append(f"{pad}{k}:")
                    emit(v, indent + 1)
                elif isinstance(v, bool):
                    lines.append(f"{pad}{k}: {'true' if v else 'false'}")
                elif isinstance(v, str):
                    lines.append(f'{pad}{k}: "{v}"')
                else:
                    lines.append(f"{pad}{k}: {v}")
        elif isinstance(obj, list):
            for item in obj:
                if isinstance(item, dict):
                    if prefix:
                        lines.append(f"{pad}{prefix}:")
                    emit(item, indent + 1)
                else:
                    lines.append(f"{pad}- {q(item)}")

    def q(v):
        if isinstance(v, str):
            return f'"{v}"'
        elif isinstance(v, bool):
            return "true" if v else "false"
        return str(v)

    emit(obj)
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# seedsGenerate: generate target integrals matching LIE's seedsGenerate[sector,r,s]
# ---------------------------------------------------------------------------
def seeds_generate(sector, r, s):
    """Generate target integral indices: corner + inddots - indrank
    sector: binary vector, r=max dot count, s=max ISP rank (negative index)"""
    n = len(sector)
    pd_indices = [i for i, v in enumerate(sector) if v == 1]  # propagator positions
    isp_indices = [i for i, v in enumerate(sector) if v == 0]  # ISP positions

    # corner = top sector base
    corner = [1 if v == 1 else 0 for v in sector]

    # inddots: all vectors of sum <= r on pd positions
    def gen_seeds(dim, level):
        res = []

        def rec(pos, remaining, current):
            if pos == dim:
                if remaining == 0:
                    res.append(list(current))
                return
            for v in range(remaining + 1):
                current[pos] = v
                rec(pos + 1, remaining - v, current)
            current[pos] = 0

        for total in range(level + 1):
            rec(0, total, [0] * dim)
        return res

    inddots = gen_seeds(len(pd_indices), r)
    indrank = gen_seeds(len(isp_indices), s)

    seeds = []
    for dot in inddots:
        for rnk in indrank:
            seed = list(corner)
            for j, idx in enumerate(pd_indices):
                seed[idx] += dot[j]
            for j, idx in enumerate(isp_indices):
                seed[idx] -= rnk[j]
            seeds.append(seed)

    # Deduplicate and sort
    seeds = [list(x) for x in set(tuple(x) for x in seeds)]
    seeds.sort(key=lambda x: (-sum(x), x))
    return seeds


# ---------------------------------------------------------------------------
# Main generation
# ---------------------------------------------------------------------------
def generate(family, output_dir, r=1, s=1):
    with open(os.path.join(FAMILIES_DIR, family + ".json")) as f:
        cfg = json.load(f)

    out = os.path.join(output_dir, family)
    config_dir = os.path.join(out, "config")
    os.makedirs(config_dir, exist_ok=True)

    prop_name = cfg.get("name", family)
    loop_mom = cfg["loopMomenta"]
    ext_mom = cfg["externalMomenta"]
    props = cfg["propagators"]
    kinematic_rules = cfg.get("kinematicRules", {})
    top_sector = cfg.get("topSector", [1] * len(props))
    numeric = cfg.get("numeric", {})
    invars = extract_invariants(kinematic_rules)

    # Determine mass symbol: only use if non-zero in numeric
    mass_symbol = None
    if "msq" in numeric and numeric["msq"] != 0:
        mass_symbol = "msq"

    # --- integralfamilies.yaml ---
    kira_props = [parse_propagator(p, mass_symbol) for p in props]
    top_sectors = [sector_to_int(top_sector)]

    family_data = {
        "integralfamilies": [
            {
                "name": prop_name,
                "loop_momenta": loop_mom,
                "top_level_sectors": top_sectors,
                "propagators": kira_props,
            }
        ]
    }
    with open(os.path.join(config_dir, "integralfamilies.yaml"), "w") as f:
        f.write(yaml_str(family_data))

    # --- kinematics.yaml ---
    scalar_rules = []
    for k, v in kinematic_rules.items():
        rule = scalar_to_kira_rule(k, v)
        if rule is not None:
            scalar_rules.append(rule)

    # mass_symbol already determined above (None if msq=0 in numeric)
    kin_invariants = [[inv, 2] for inv in invars if inv != "msq"]
    if mass_symbol and mass_symbol not in [x[0] for x in kin_invariants]:
        kin_invariants.append([mass_symbol, 2])

    replace_one = None
    for inv in invars:
        if inv == "msq":
            continue
        num_val = numeric.get(inv, None)
        if num_val is not None and num_val != 0:
            replace_one = inv
            break
    if replace_one is None and mass_symbol:
        replace_one = mass_symbol
    if replace_one is None and kin_invariants:
        replace_one = kin_invariants[0][0]

    kinematics_data = {
        "kinematics": {
            "incoming_momenta": ext_mom,
            "outgoing_momenta": [],
            "momentum_conservation": [],
            "kinematic_invariants": kin_invariants,
            "scalarproduct_rules": scalar_rules,
        }
    }
    if replace_one:
        kinematics_data["kinematics"]["symbol_to_replace_by_one"] = replace_one

    with open(os.path.join(config_dir, "kinematics.yaml"), "w") as f:
        f.write(yaml_str(kinematics_data))

    # --- jobs.yaml ---
    jobs_data = {
        "jobs": [
            {
                "reduce_sectors": {
                    "reduce": [{"sectors": top_sectors, "r": r, "s": s}],
                    "select_integrals": {
                        "select_mandatory_list": [[prop_name, "integrals"]]
                    },
                    "run_initiate": True,
                    "run_triangular": True,
                    "run_back_substitution": True,
                }
            },
            {"kira2math": {"target": [[prop_name, "integrals"]]}},
        ]
    }
    with open(os.path.join(out, "jobs.yaml"), "w") as f:
        f.write(yaml_str(jobs_data))

    # --- integrals file ---
    # Generate target integrals matching LIE's seedsGenerate[sector, r, s]
    seeds = seeds_generate(top_sector, r, s)
    integrals = [f"{prop_name}[{','.join(str(x) for x in s)}]" for s in seeds]

    with open(os.path.join(out, "integrals"), "w") as f:
        f.write("\n".join(integrals) + "\n")

    print(f"  Generated: {out}/")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    import argparse
    parser = argparse.ArgumentParser(description="Generate Kira input from LIE family JSON")
    parser.add_argument("families", nargs="*", help="Family names")
    parser.add_argument("--all-1l2p", action="store_true", help="Generate all 1L + 2L2P families")
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT, help="Output directory")
    parser.add_argument("--r", type=int, default=1, help="Maximum dot count on propagators (default: 1)")
    parser.add_argument("--s", type=int, default=1, help="Maximum ISP rank / negative index count (default: 1)")
    args = parser.parse_args()

    targets = []
    if args.all_1l2p:
        targets = ALL_1L2P
    elif args.families:
        targets = args.families
    else:
        parser.print_help()
        sys.exit(1)

    for fam in targets:
        generate(fam, args.output_dir, r=args.r, s=args.s)

    print(f"\nDone. {len(targets)} families generated in {args.output_dir}/ (r={args.r}, s={args.s})")


if __name__ == "__main__":
    main()
