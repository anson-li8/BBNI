# BBNI Function Analysis

---
## `update.ancestor_matrix`

**What it does:** It calculates the transitive closure of the network's directed adjacency matrix to locate all ancestor-descendant relationships. This helps in ensuring the proposed MCMC topologies represent directed acyclic graphs (DAG) without cyclic loops.

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

**What it does:** It ensures that the proposed network topology fulfills the DAG constraint by checking the diagonal (self-referring) of ancestor matrix for illegal cyclic loops.

**Arguments:** - `ances_matrix` — (matrix) The transitive closure matrix of the network topology 𝑇, where (i, j) = 1 means node i is an ancestor of node j

**Returns:** (numeric) Integer count of the number of nodes that are ancestors of themselves. If greater than 0, the proposed topology has cyclic loops and must be rejected.

**Hidden globals:** None

**Paper reference:** "The MCMC Algorithm" section, specifically where it talks about the constraints on network topology updates: "There are three types of MCMC moves to update the parent sets: adding parent(s), removing parent(s) and swapping parent(s)... But if adding parent(s) leads to a cyclic graph, that specific move is illegal". This function helps enforce that acyclic requirement by returning a loop count greater than 0 to indicate a topology that must be rejected.

**Status:** [x] Analyzed [ ] Cleaned [ ] Documented

---
## `Prop_Trans_Func_Matrix`

**What it does:** For a proposed network topology, assigns an integer code to each nonzero entry that indicates the Boolean function used by the child node. The code range is based on parent count (1 or 2): Codes 1–10 correspond to ten non-degenerate two-input Boolean functions, while codes 11–12 correspond to the unary functions of identity and negation.

**Arguments:**
- `prop_incid_matrix` — (matrix) Proposed incidence matrix for network topology

**Returns:**
Transition-function matrix with same dimensions as input. 0 indicates no incoming edge. Non-zero entries each contain a code from 1-12 indentifying which Boolean operation each genes utilizes for calculation of the input of its parents(s).

**Hidden globals:** 
None

**Paper reference:** "Prior Distributions" subsection, specifically where it explains the allowable Boolean update rules: "If $W(g_i)$ is the set {a}, $f_i(a)$ can be either a or $\overline{a}$; if $W(g_i)$ is the set {a, b}, $f_i(a,b)$ has 10 non-degenerative choices...". This function follows that logic by mapping each parent set to its corresponding valid Boolean transition rule.

**Status:** [x] Analyzed [ ] Cleaned [ ] Documented

---
## `Error_LLH`
**What it does:** It calculates the collapsed posterior log-probability of a proposed network topology (T) and its Boolean functions (F) given the observed data (G). It parses the TRFUM into root and non-root nodes, predicts each node's value based on its Boolean rule, counts mismatches, and evaluates the Beta-Binomial collapsed likelihood. By integrating out Bernoulli parameters $θ=\{𝑝_𝐸,𝑝_𝑖\}$ under Beta priors, the function returns $logp(T,F|G)$, which is the model score for the MCMC algorithm.

**Arguments:**
- `TRFUM` - (matrix) Transition Function Matrix combining both the current network topology T and the specific Boolean logic functions F assigned to each node (integer codes 1-12).

**Returns:**(list) List of two vectors containing evaluating model fit. results[[1]] contains [ErrorFactor, RootFactor, likelihood, post_para, log_post_model], with the most notable being log_post_model, which is the collapsed posterior metric $logp(T,F|G)$. results[[2]] contains [para_sample, mismath, Perror], containing point estimates of the root ON-probabilities, total mismatches, and estimate noise error rate, e.

**Hidden globals:**
- `GeneData` - (matrix) Observed binary gene expression time-series dataset G (rows = nodes, columns = time points) 
- `SampleSize` - (integer) Total number of time points/samples in GeneData 
- `num.node` - (integer) Total number of genes/nodes in the network 
- `prior_para` - (matrix) Beta distribution hyperparameters ($\alpha, \beta$) for root nodes and the global noise parameter $e$ 
- `penalty` - (numeric) The structural prior probability per edge used to penalize network complexity $P(T)$

**Paper reference:**“Posterior Distributions” subsection: “To circumvent this problem, we analytically integrate out all pᵢ’s and $p_E$ from the above posterior distribution, which results in the following collapsed version of the posterior distribution: Equation (4)” “We have designed an MCMC algorithm to sample from $𝑝(𝐹,𝑇∣𝐺)$, which avoids the dimension change caused by pᵢ’s.” This function directly implements the collapsed posterior model described in Equation (4), using mismatch counts 𝐵 and root ON‑counts $C_i$ to compute the integrated likelihood.

**Status:** [x] Analyzed [ ] Cleaned [ ] Documented

---
## `GenerateNetwork`
**What it does:** Randomly generates an initial, legal DAG graph topology $T$ and assigns a corresponding Boolean transition function $F$ to each node, ensuring a maximum in-degree value of 2.

**Arguments:**
- `num.node` - (integer) Total number of genes/nodes in the network

**Returns:**(matrix) Square transition function matrix combining both the initial DAG topology and the randomly assigned Boolean logic functions (integer code from 1-12). Ensures graph is acyclic before returning

**Hidden globals:**
- None. Relies on `update.ancestor_matrix` and `check.ances.matrix` functions to verify DAG constraint

**Paper reference:**"Model" subsection, stating structural constraints of the Boolean network: “we will focus on the case where the maximum in-degree of all nodes in the network is bounded by 2” and that the required structure is a “directed acyclic graph denoted by ${G, T, F}$.”  This function generates an initial topology that fulfills those requirements by assigning random parent sets of size 0-2 and rejecting any topology that forms a directed loop.

**Status:** [x] Analyzed [ ] Cleaned [ ] Documented

---
## `GenerateSample`
**What it does:** Simulate a time-series observation dataset $G$ by identifying root nodes to simulate independently, then topologically computing the remaining nodes from $t-1$ to time $t$ using their specified Boolean logic functions and an inverted noise parameter $e$.

**Arguments:**
- `trans_matrix` - (matrix) Square matrix combining network topology $T$ and integer coded Boolean logic functions $F$ assigned to each directed edge.

**Returns:**(matrix) Simulated binary gene expression matrix $G$, where rows represent individual genes/nodes and columns represents sequential points in time, acting as final output data with simulated biological noise

**Hidden globals:**
- `num.node` - (integer) Total number of network nodes
- `SampleSize` - (integer) Total number of time points to simulate
- `para` - (numeric vector) Baseline success probabilities $\theta_i$ used to generate the expression states of root nodes by independent Bernoullli trials.
- `error` - (matrix) Pre-generated binary noise matrix. Applied by `bitXor` operation to occasionally flip Boolean output, injecting natural noise expected by model

**Paper reference:**
"Model" subsection, where the data‑generation equations are defined. Root nodes follow the Bernoulli model in Equation (1):
“If $W(g_i)$ is an empty set... we assume an independent Bernoulli distribution for it, i.e., 
$Pr(g_{ij}=1)=p_i$.” Non‑root nodes follow the noisy Boolean update in Equation (1): “$g_i=f_i(W(g_i))\oplus\epsilon$”, where $\epsilon$ is a Bernoulli noise term with probability $p_\epsilon$ of flipping the output. This function implements these equations by randomizing roots independently and applying corresponding Boolean logic with XOR‑calculated noise for the non‑roots.

**Status:** [x] Analyzed [ ] Cleaned [ ] Documented

---
## `BF1` through `BF14`
**What it does (combined):** Each function (BF1-BF10) for two-input functions, BF11-BF14 for single-input) computes a Bayesian Information Criterion (BIC) score for a candidate Boolean function by counting state mismatches between predicted and observed outputs. Lower BIC means better fit. The function calculates this efficiently by compressing the raw data into compact 3-gene frequency count (truth table) and using it to instantly find matching patterns.

**Arguments:**
- `test.stat` - (integer vector) Frequency count vector [c000, c001, c010, c011, c100, c101, c110, c111] storing how many times each specific combination of parent and child states occurs together in the data

**Returns:**
Numeric BIC value (or NA if mismatch count exceeds threshold). Used to rank candidate functions.

**Hidden globals:**
- `SampleSize` - (integer) Sample size
- `pseudo.count` - (numeric) Pseudocount to avoid zero cells
- `threshold` - (numeric) Rejection threshold for models with high error rates

**Paper reference:**
"Prior Distributions" subsection, specifically where it explains the allowable Boolean update rules: "If $W(g_i)$ is the set {a}, $f_i(a)$ can be either a or $\overline{a}$; if $W(g_i)$ is the set {a, b}, $f_i(a,b)$ has 10 non-degenerative choices...". These 14 functions (10 two-input + 4 one-input) are the candidates evaluted by the BIC scorign functions BF1-BF14.

**Function mappings:**
- `BF1-BF10`: Two-input Boolean functions (AND, NAND, OR, NOR, OR-NOT, NOT-OR, AND-NOT, NOT-AND, XOR, NXOR)
- `BF11-BF12`: Pairwise relations (identity, identity of second parent)
- `BF13-BF14`: Negations (NOT of first, NOT of second)

**Status:** [x] Analyzed [ ] Cleaned [ ] Documented

---
## `ProposalConstruction`
**What it does:** Pre-computes a weighted proposal distribution for the MCMC algorithm. It evaluates all possible 1-input and 2-input Boolean logic functions using BIC. The inverse BIC of the potential candidates are used as selection weights for the MCMC sampler.

**Arguments:**
- `GeneData` - (matrix) Observational gene expression data

**Returns:** (list) Containing two matrices (Candidate[[1]] for 2-input triplets and Candidate[[2]] for 1-input pairs). These matrices store the candidate parent nodes, the best-fitting Boolean logic function, the raw miscounts, and the calculated proposal selection weights.

**Hidden globals:**
- SampleSize — (integer) Total number of time points. 
**CRITICAL NOTE:** There is already a defined local variable sample.size <- ncol(gene.data) at the top of the function, but the global variable SampleSize is used instead inside the nested loop during subsetting: gene.data[i,1:(SampleSize-1)]. 
- BF1 through BF14 — (functions) Fourteen external helper functions called in the inner loop to evaluate the BICS of the transition state counts (test.stat) and return the total data mismatches for each specific Boolean rule. (Note: Rules 11-14 in this script map to the 1-input rules for $g_i$, $g_j$, and their complements).

**Paper reference:**
"Prior Distributions" subsection, specifically where it explains the allowable Boolean update rules: "If $W(g_i)$ is the set {a}, $f_i(a)$ can be either a or $\overline{a}$; if $W(g_i)$ is the set {a, b}, $f_i(a,b)$ has 10 non-degenerative choices...". These 14 functions are the candidates being evaluated and weighted in this distribution construction step.

**Status:** [x] Analyzed [ ] Cleaned [ ] Documented