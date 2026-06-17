#!/usr/bin/env python3
"""Compare MMA and C++ LIE relation subspaces."""
import re, sys

MODULUS = 179424673

# ============ C++ Parser ============
def extract_braced(text, idx):
    i = text.find('{', idx)
    if i < 0: return None
    d = 0
    for j in range(i, len(text)):
        if text[j] == '{': d += 1
        elif text[j] == '}':
            d -= 1
            if d == 0: return text[i+1:j]
    return None

def parse_cpp(fp):
    with open(fp) as f:
        txt = f.read()
    for key in ['Alphas', 'Betas', 'Coefficients']:
        i = txt.find(f'"{key}"'); a = txt.find('->', i)
        if key == 'Alphas':
            alphas = [(int(m.group(1)),int(m.group(2))) for m in re.finditer(r'\{(\d+),(\d+)\}', extract_braced(txt, a))]
        elif key == 'Betas':
            betas = [(int(m.group(1)),int(m.group(2))) for m in re.finditer(r'\{(\d+),(\d+)\}', extract_braced(txt, a))]
        else:
            block = extract_braced(txt, a)
            rows = [[int(x) for x in re.findall(r'\d+', row)] for row in block.split('\n') if '{' in row]
    # Skip first column (Particular solution)
    ab = [(a, b) for a in alphas for b in betas]
    sols = []
    for si in range(1, len(rows[0])):
        vec = {}
        for vi, row in enumerate(rows):
            v = row[si]
            if v != 0: vec[ab[vi]] = v % MODULUS
        sols.append(vec)
    return sols

# ============ MMA Parser ============
def parse_mma_term(term):
    """Parse a single MMA expression term."""
    term = term.strip()
    if not term: return None
    sign = -1 if term.startswith('-') else 1
    term = term.lstrip('+-').strip()
    gm = re.search(r'g\[([^\]]+)\]', term)
    if not gm: return None
    args = [a.strip() for a in gm.group(1).split(',')]
    alphas = []
    for a in args:
        m = re.search(r'v\d\s*-\s*(\d+)', a)
        if m: alphas.append(int(m.group(1))); continue
        m = re.search(r'[-]?(\d+)\s*[+-]\s*v\d', a)
        if m:
            if a.lstrip().startswith('-'): alphas.append(int(m.group(1)))
            else: alphas.append(-int(m.group(1)))
            continue
        alphas.append(0)
    before = term[:gm.start()].strip()
    if not before:
        return (tuple(alphas), (0,0), sign)
    coeff = sign; b1 = b2 = 0
    for p in re.split(r'\*', before):
        p = p.strip()
        mm = re.match(r'v1\^(\d+)$', p)
        if mm: b1 += int(mm.group(1)); continue
        mm = re.match(r'v2\^(\d+)$', p)
        if mm: b2 += int(mm.group(1)); continue
        if p == 'v1': b1 += 1; continue
        if p == 'v2': b2 += 1; continue
        if p in ('v1*v2','v2*v1'): b1 += 1; b2 += 1; continue
        try: coeff *= int(p)
        except: pass
    return (tuple(alphas), (b1,b2), coeff % MODULUS)

def parse_mma_rel(expr):
    """Parse MMA relation expression - bracket-aware split by + and -."""
    s = expr.strip('{}').strip()
    # Split by + or - at bracket-level 0
    terms = []
    depth = 0
    cur = []
    prev_was_plus = True  # If first char is -, it's a sign, not an operator
    for c in s:
        if c in '[(': 
            depth += 1
        elif c in '])': 
            depth -= 1
        elif (c == '+' or c == '-') and depth == 0:
            t = ''.join(cur).strip().lstrip('+').strip()
            if t: terms.append(t)
            cur = []
            if c == '-': 
                cur.append('-')  # Keep the minus sign for the next term
            continue
        cur.append(c)
    if cur:
        t = ''.join(cur).strip().lstrip('+').strip()
        if t: terms.append(t)
    result = {}
    for t in terms:
        r = parse_mma_term(t)
        if r:
            key = (r[0], r[1])
            result[key] = (result.get(key,0) + r[2]) % MODULUS
    return result

def load_mma(fp):
    with open(fp) as f:
        txt = f.read()
    txt = txt.replace('"', '')
    m = re.search(r'\$MMARelations\s*=\s*(\{.*)', txt, re.DOTALL)
    if not m: return []
    body = m.group(1).rstrip(';').strip()
    # Find innermost {groups} with g[ that have NO nested braces
    # {[^{}]*g\[[^{}]*} matches {content} where content has g[ but no nested braces
    rels = []
    for m in re.finditer(r'\{([^{}]*g\[[^{}]*)\}', body):
        inner = m.group(1).strip()
        if not inner: continue
        # Split by commas at bracket-level 0 (handling g[...] brackets)
        parts = []; dep = 0; cur = []
        for c in inner:
            if c in '[(': dep += 1
            elif c in '])': dep -= 1
            elif c == ',' and dep == 0:
                parts.append(''.join(cur).strip()); cur = []; continue
            cur.append(c)
        if cur: parts.append(''.join(cur).strip())
        for p in parts:
            p = p.strip().strip(',').strip()
            if p and 'g[' in p:
                rels.append(p)
    return rels

# ============ Subspace Comparison ============
def rank_mod(rows, p):
    if not rows or not rows[0]: return 0
    m = [[x%p for x in r] for r in rows]
    n, c = len(m), len(m[0]); rk = 0
    for col in range(c):
        pi = next((i for i in range(rk, n) if m[i][col]%p!=0), None)
        if pi is None: continue
        m[rk], m[pi] = m[pi], m[rk]
        inv = pow(m[rk][col], -1, p)
        for k in range(col, c): m[rk][k] = (m[rk][k]*inv) % p
        for i in range(n):
            if i != rk and m[i][col] != 0:
                f = m[i][col]
                for k in range(col, c): m[i][k] = (m[i][k]-f*m[rk][k]) % p
        rk += 1
    return rk

if __name__ == '__main__':
    base = "/root/Large-Index-Expansion-CPP/verify/bub00"
    
    # C++
    cpp_all = []
    for deg in [0, 1, 2]:
        sols = parse_cpp(f"{base}/Compare-CPPRelation-bub00_lev2_deg{deg}.m")
        cpp_all.extend(sols)
        print(f"C++ deg={deg}: {len(sols)} basis vectors")
    print(f"C++ total: {len(cpp_all)}")
    
    # MMA
    raw_rels = load_mma(f"{base}/Compare-MMARelation-lev2.m")
    mma_all = [parse_mma_rel(r) for r in raw_rels]
    mma_all = [d for d in mma_all if d]
    print(f"MMA total: {len(mma_all)}")
    
    if not cpp_all or not mma_all:
        print("ERROR: No data!"); sys.exit(1)
    
    # Common basis
    keys = set()
    for d in mma_all: keys.update(d.keys())
    for d in cpp_all: keys.update(d.keys())
    def sk(k): a,b=k; return (sum(a),a[0],a[1],sum(b),b[0],b[1])
    keys = sorted(keys, key=sk)
    print(f"\nBasis size: {len(keys)}")
    for k in keys[:6]: print(f"  α={k[0]} β={k[1]}")
    if len(keys) > 6: print(f"  ... + {len(keys)-6} more")
    
    mma_mat = [[d.get(k,0)%MODULUS for k in keys] for d in mma_all]
    cpp_mat = [[d.get(k,0)%MODULUS for k in keys] for d in cpp_all]
    
    r1 = rank_mod(mma_mat, MODULUS)
    r2 = rank_mod(cpp_mat, MODULUS)
    r3 = rank_mod(mma_mat + cpp_mat, MODULUS)
    
    print(f"\n  rank(MMA)  = {r1}")
    print(f"  rank(C++)  = {r2}")
    print(f"  rank(comb) = {r3}")
    
    if r1 == r2 == r3:
        print("\n✅ 子空间一致！MMA 和 C++ 产生相同的线性关系子空间")
    elif r3 > max(r1, r2):
        print(f"\n❌ 子空间不一致：维度差 = {r3 - max(r1, r2)}")
    else:
        print("\n⚠️ 部分一致")
