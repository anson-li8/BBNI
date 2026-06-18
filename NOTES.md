# BBNI Function Analysis
---
## `update.ancestor_matrix`

**What it does:** It calculates the transitive closure of the network's directed adjacency matrix to locate all ancestor-descendant relationships. This helps in ensuring the proposed MCMC topologies represent directed acyclic graphs without cyclic loops.

**Arguments:**
- `incid_matrix` — (matrix) A square binary adjacency matrix that represents the current network topology 𝑇

**Returns:** 
A binary matrix of the same dimensions as incid_matrix. Entry (i, k) = 1 means that node 1 is an ancestor of node k through one or more directed edges.

**Hidden globals:** 
None

**Paper reference:** "The MCMC Algorithm" section, specifically where it talks about the constraints on network topology updates: "There are three types of MCMC moves to update the parent sets: adding parent(s), removing parent(s) and swapping parent(s)... But if adding parent(s) leads to a cyclic graph, that specific move is illegal". This function outputs to the algorithm which nodes are upstream of which others, so it can prevent cycles when updating the network.

**Status:** [x] Analyzed [ ] Cleaned [ ] Documented

---
## `check.ances.matrix`

**What it does:** It ensures that the proposed network topology fulfills the directed acyclic graph constraint by checking the diagonal (self-referring) of ancestor matrix for illegal cyclic loops.

**Arguments:**
- `ances_matrix` — (matrix) The transitive closure matrix of the network topology 𝑇, where (i, j) = 1 means node i is an ancestor of node j

**Returns:** (numeric)
Integer count of the number of nodes that are ancestors of themselves. If greater than 0, the proposed topology has cyclic loops and must be rejected.

**Hidden globals:** 
None

**Paper reference:** "The MCMC Algorithm" section, specifically where it talks about the constraints on network topology updates: "There are three types of MCMC moves to update the parent sets: adding parent(s), removing parent(s) and swapping parent(s)... But if adding parent(s) leads to a cyclic graph, that specific move is illegal". This function helps enforce that acyclic requirement by returning a loop count greater than 0 to indicate a topology that must be rejected.

**Status:** [x] Analyzed [ ] Cleaned [ ] Documented

---
## `Prop_Trans_Func_Matrix`

**What it does:** For a proposed network topology, assigns an integer code to each nonzero entry that indicates the Boolean function used by the child node. The code range is based on parent count (1 or 2): Codes 1–10 correspond to ten non-degenerate two-input Boolean functions, while codes 11–12 correspond to the unary functions of identity and negation.

**Arguments:**
- `prop_incid_matrix` — (matrix) Proposed incidence matrix for network topology

**Returns:**
Matrix with same dimensions as input. Zeros indicate no incoming edge. Non-zero entries each contain a code from 1-12 indentifying which Boolean operation each genes utilizes for calculation of the input of its parents(s).

**Hidden globals:** 
None

**Paper reference:** "Prior Distributions" subsection, specifically where it explains the allowable Boolean update rules: "If W(g_i) is the set {a}, f_i(a) can be either a or ¬a; if W(g_i) is the set {a, b}, f_i(a,b) has 10 non-degenerative choices...". This function follows that logic by mapping each parent set to its corresponding valid Boolean transition rule.

**Status:** [x] Analyzed [ ] Cleaned [ ] Documented

---
## `Error_LLH`

**What it does:** 

**Arguments:**
- 

**Returns:**

**Hidden globals:** 

**Paper reference:** 

**Status:** [x] Analyzed [ ] Cleaned [ ] Documented