#' Calculate Transitive Closure of Adjacency Matrix
#'
#' Calculates the transitive closure of the network's directed adjacency matrix
#' to locate all ancestor-descendant relationships. This is used by the MCMC
#' algorithm to ensure proposed network topologies represent directed acyclic
#' graphs (DAGs) and to prevent illegal cyclic moves when adding or swapping parents.
#' @param incid_matrix A square binary adjacency matrix representing the current network topology.
#'
#' @returns A binary matrix of the same dimensions as `incid_matrix`. An entry of 1 at (i, k) indicates that node i is an ancestor of node k through one or more directed edges.
#' @noRd
update_ancestor_matrix <- function(incid_matrix) #  for a given incidence matrix, generate corresponding ancestor matrix
{
  ances_matrix <- matrix(0, nrow = nrow(incid_matrix), ncol = ncol(incid_matrix)) # important! every time before updating, ancestor matrix should be cleared.
  ances_matrix <- incid_matrix
  for (ii in 1:nrow(ances_matrix))
  {
    for (i in 1:nrow(ances_matrix)) {
      for (j in 1:ncol(ances_matrix)) {
        for (k in 1:nrow(ances_matrix)) {
          if (ances_matrix[i, j] == 1 && ances_matrix[j, k] == 1) { # it should be updated num.node times
            ances_matrix[i, k] <- 1
          }
        }
      }
    }
  }
  return(ances_matrix)
}

#' Check Ancestor Matrix for Cyclic Loops
#'
#' Ensures that the proposed network topology fulfills the directed acyclic graph
#' (DAG) constraint by checking the diagonal of the ancestor matrix for illegal
#' cyclic loops. This enforces the requirement that adding or swapping parents
#' cannot lead to a cyclic graph during MCMC topology updates.
#'
#' @param ances_matrix A binary transitive closure matrix of the network topology, where an entry of 1 at (i, j) means node i is an ancestor of node j.
#'
#' @return An integer count of the number of nodes that are ancestors of themselves. A value greater than 0 indicates the proposed topology has cyclic loops and must be rejected.
#' @noRd
check_ances_matrix <- function(ances_matrix) # check whether there are loops in the whole network by checking ancestor matrix
{
  loop <- 0
  for (i in 1:nrow(ances_matrix)) {
    if (ances_matrix[i, i] == 1) {
      loop <- loop + 1
    }
  }
  return(loop)
}
