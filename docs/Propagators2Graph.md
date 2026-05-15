# From Propagators (Edge-Cycle Matrix) to Graph Topology

**Problem**: Given only the edge-cycle matrix $C \in \{0, \pm1\}^{c \times m}$ ($c = m - n + 1$ for a connected graph) — equivalently, the set of directed cycles in terms of which edges they traverse and with what relative sign — reconstruct the graph topology and, in particular, a vertex-star basis of the cut space. This is the algebraic inverse of constructing a cycle matrix from a graph.

---

## 1. Algorithm Overview

### Step 1: Column reduction → spanning tree

Perform **Gaussian elimination on columns** of $C$:
- Pivot on columns of $C$ to obtain a column-echelon form
- The $c$ pivot columns are the **chords** (co-tree edges) $R$, $|R| = c$
- The remaining $n-1$ non-pivot columns are the **tree edges** $T$, $|T| = n-1$

From the echelon form we also read off the coefficients $\alpha_{fe} \in \{0, \pm1\}$ expressing each chord column in terms of tree-edge columns:

$$ C_f = \sum_{e \in T} \alpha_{fe} \, C_e, \qquad f \in R $$

The **fundamental cycle** of chord $f$ is:
$$ C_f = \{\,e \in T : \alpha_{fe} \neq 0\,\} \cup \{f\} $$

### Step 2: Fundamental cuts

For each tree edge $e \in T$, the **fundamental cut** $\delta(T_e)$ is:
$$
\delta(T_e)_g =
\begin{cases}
+1 \text{ (or } -1) & g = e \\
\alpha_{fe} & g = f \in R,\ \alpha_{fe} \neq 0 \\
0 & \text{otherwise}
\end{cases}
$$

The set $\{\delta(T_e) : e \in T\}$ is a $\{0,\pm1\}$-basis of the cut space.

### Step 3: Tree reconstruction from fundamental cuts

The key observation: two tree edges $e_1, e_2 \in T$ share a vertex in the tree $T$ **if and only if** a specific condition on their chord-support sets holds (see §2 below). Using this, we reconstruct the tree $T$ combinatorially.

### Step 4: Vertex-star basis

Once the tree topology is known, pick a root vertex $v_0$. For each vertex $v \neq v_0$, let $\text{path}_T(v_0, v)$ be the unique path in $T$. Then:

$$ \delta(v) = \sum_{e \, \in \, \text{path}_T(v_0, v)} \delta(T_e) $$

The summation is ordinary integer addition (no cancellation since the edges are distinct). The set $\{\delta(v) : v \in V \setminus \{v_0\}\}$ is the **vertex-star basis**.

---

## 2. Tree-Edge Adjacency Criterion

Let $S(e) := \{f \in R : \alpha_{fe} \neq 0\}$ be the set of chords whose fundamental cycle includes tree edge $e$.

**Theorem** (Tree-edge adjacency). Let $G$ be 2-connected and simple. Two tree edges $e_1, e_2 \in T$ share a vertex in $T$ **iff** both:
1. $S(e_1) \cap S(e_2) \neq \varnothing$
2. $\nexists\, e_3 \neq e_1, e_2$ such that $S(e_3) \supseteq S(e_1) \cap S(e_2)$

**Proof.**

Let endpoints of chord $f$ be $p_f, q_f$. Then $f \in S(e)$ iff $p_f$ and $q_f$ lie in different components of $T - e$.

*($\Rightarrow$)* Let $e_1 = (a,b)$, $e_2 = (b,c)$ share vertex $b$. Since $G$ is 2-connected, there exists a chord $f$ connecting a vertex in the component of $T - e_1$ that does not contain $b$ to a vertex in the component of $T - e_2$ that does not contain $b$. This chord necessarily belongs to $S(e_1) \cap S(e_2)$, verifying condition (1).

For condition (2): any $e_3 \neq e_1, e_2$ either lies off the path between the far ends of $e_1, e_2$ (in which case a chord exists that avoids $e_3$), or lies on the path — but the path only contains $e_1, e_2$, so no such $e_3$ exists.

*($\Leftarrow$)* Suppose $e_1, e_2$ do not share a vertex.

**Case A:** $e_1, e_2$ lie in different branches from a branching vertex. Then any $T$-path of a chord can contain at most one of them, so $S(e_1) \cap S(e_2) = \varnothing$, violating (1).

**Case B:** $e_1, e_2$ are collinear in $T$ with at least one intermediate edge $e_3$ between them, i.e., $e_1 - e_3 - \dots - e_k - e_2$, $k \ge 1$. By 2-connectivity, $S(e_1) \cap S(e_2) \neq \varnothing$, so (1) holds. But any chord that crosses both $e_1$ and $e_2$ must traverse the entire path between them, hence must cross every intermediate edge, giving $S(e_1) \cap S(e_2) \subseteq S(e_3)$. This violates (2).

Thus both conditions hold iff $e_1, e_2$ are adjacent in $T$. $\square$

---

### 2.1 Line graph of $T$ from chord intersections

Define the **tree-edge intersection matrix**:
$$ M_{e_1 e_2} = |S(e_1) \cap S(e_2)|, \qquad e_1, e_2 \in T $$

The adjacency graph of $T$ (its line graph) is recovered from $M$: two edges are adjacent in $T$ iff $M_{e_1 e_2} > 0$ and there is no $e_3$ with $M_{e_1 e_2} \le M_{e_1 e_3}$ and $M_{e_1 e_2} \le M_{e_3 e_2}$ satisfying the strict inclusion of support sets above.

A tree is uniquely determined by its line graph (unless $T \cong K_{1,3}$, whose line graph is $K_3$ and could also arise as the line graph of $K_3$ — but $K_3$ is not a tree). This edge case can be resolved by checking the intersection pattern directly.

---

## 3. Tree Reconstruction Algorithm

```
Input: fundamental cuts δ(T_e), chord coefficients α_{fe}
Output: tree T (vertex-labeled), path matrix D

1. For each e ∈ T, compute S(e) = {f ∈ R : α_{fe} ≠ 0}

2. Compute adjacency graph L(T) on T:
   For each pair e₁, e₂ ∈ T:
       if S(e₁) ∩ S(e₂) ≠ ∅ AND
          ¬∃ e₃ ≠ e₁, e₂ with S(e₃) ⊇ S(e₁) ∩ S(e₂):
           mark e₁—e₂ as adjacent in L(T)

3. Reconstruct T from its line graph L(T):
   a. Each vertex of L(T) is an edge of T
   b. Each edge of L(T) means two tree edges share a vertex
   c. Reconstruct the incidence of T by taking L(T), picking a
      maximal clique for a vertex, and growing outward
   d. Since T is a tree, its line graph L(T) is chordal and the
      reconstruction is unique (with the K₁,₃ caveat above)

4. Pick root vertex v₀ on T
5. For each vertex v, let D_{v,e} = 1 iff e lies on the unique
   path from v₀ to v in T
6. Compute vertex-star basis:
       δ(v) = Σ_e D_{v,e} · δ(T_e)
```

---

## 4. Algebraic Shortcut (Matrix Form)

Arrange fundamental cuts as rows of matrix $K \in \{0,\pm1\}^{(n-1)\times m}$.

Let $B_T \in \{0,\pm1\}^{(n-1)\times(n-1)}$ be the oriented incidence matrix of tree $T$ (rows = non-root vertices, columns = tree edges). $B_T$ is invertible and its inverse has the combinatorial interpretation:

$$ (B_T^{-1})_{v,e} = \begin{cases}
+1 & \text{if } e \text{ lies on the oriented path } v_0 \to v \text{ with same orientation} \\
-1 & \text{if } e \text{ lies on the path with opposite orientation} \\
0  & \text{otherwise}
\end{cases} $$

Then the vertex-star matrix $S \in \{0,\pm1\}^{n \times m}$ (with root row omitted) is:

$$ S = B_T^{-1} \cdot K $$

So the vertex-star basis is just a linear transform (the path-sign matrix) applied to the fundamental cuts.

---

## 5. Uniqueness Conditions

The vertex-star basis constructed above is unique **up to choice of root vertex** $v_0$ (i.e., which row is omitted) **iff**:

| Condition | Consequence |
|-----------|-------------|
| Graph is **simple** (no multi-edges, no self-loops) | Supports of $\delta(u)$ and $\delta(v)$ intersect in at most $\{uv\}$ |
| Graph is **2-connected** | Every tree edge has $S(e) \neq \varnothing$; bridges cannot be oriented |
| Graph is **3-connected** | The graph is uniquely determined from $C$ up to isomorphism (Whitney 2-isomorphism theorem) |

For non-3-connected graphs, $C$ determines the graph only up to 2-isomorphism, and multiple vertex-star bases exist.

---

## 6. Relation to Syzmanszik Polynomial

The Syzmanszik polynomial $S(t)$ encodes the distribution of edge weights in fundamental cycles. Specifically, for tree edge $e$, the number of chords in $S(e)$ is the degree of $e$ in the Syzmanszik polynomial sense. The inclusion poset $\{S(e) : e \in T\}$ partially ordered by subset inclusion is isomorphic to the tree $T$ itself, which provides an alternative route to reconstruction: build the Hasse diagram of the $S(e)$ poset, then read off the tree adjacency.

---

## References

- Whitney, H. "2-Isomorphic Graphs." *Amer. J. Math.* 55 (1933), 245–254.
- Syzmanszik, private communication on propagator-cycle duality in Feynman integral reductions.
