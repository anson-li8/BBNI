# BBNI Function Analysis
---
## `update.ancestor_matrix`
**Line:** 7
**File:** bbni.R

**What it does:** It calculates the transitive closure of the network's directed adjacency matrix to locate all ancestor-descendant relationships. This ensures the proposed MCMC topologies represent directed acyclic graphs without cyclic loops.

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
**Line:** 22
**File:** bbni.R

**What it does:** It ensures that the proposed network topology fulfills the directed acyclic graph constraint by checking the diagonal (self-referring) of ancestor matrix for illegal cyclic loops.

**Arguments:**
- `ances_matrix` — (matrix) The transitive closure matrix of the network topology 𝑇, where (i, j) = 1 means node i is an ancestor of node j

**Returns:** (numeric)
Integer count of the number of nodes that are ancestors of themselves. If greater than 0, the proposed topology has cyclic loops and must be rejected.

**Hidden globals:** 
None

**Paper reference:** "The MCMC Algorithm" section, specifically where it talks about the constraints on network topology updates: "There are three types of MCMC moves to update the parent sets: adding parent(s), removing parent(s) and swapping parent(s)... But if adding parent(s) leads to a cyclic graph, that specific move is illegal". This function helps enforce that acyclic requirement by returning a loop count greater than 0 to indicate a topology that must be rejected.

**Status:** [x] Analyzed [ ] Cleaned [ ] Documented