#' Execute Metropolis-within-Gibbs MCMC Sampler for Boolean Networks
#'
#' Executes a Metropolis-within-Gibbs Markov chain Monte Carlo (MCMC) algorithm
#' to sample from the joint posterior distribution of network topologies ($T$)
#' and Boolean logic transition functions ($F$). The sampler loops through
#' individual network nodes and proposes parent set mutations (additions,
#' removals, or swaps) among 14 candidate Boolean rules. Proposed states are
#' verified to follow the directed acyclic graph (DAG) constraint and evaluated
#' with a Metropolis-Hastings acceptance gate using the log-posterior ratios.
#'
#' @param GeneData A matrix of the observational binary expression data ($G$).
#' @param num.node An integer representing the total number of network nodes.
#' @param SampleSize An integer representing the total number of time points in the dataset.
#' @param prior_para A matrix of Beta prior hyperparameters \eqn{\alpha} and \eqn{\beta} for root node probabilities and the global noise parameter e.
#' @param num_update An integer representing the total number of MCMC iterations to perform.
#' @param penalty A numeric value representing the structural prior probability per edge used to penalize network complexity $P(T)$.
#' @param prop.ratio A numeric probability threshold used to decide whether to sample a move from the empirical proposal distribution or a uniform random distribution.
#'
#' @return A list containing the full trajectory of the MCMC chain. Specifically, `networks` (a list of sampled transition function matrices) and `log_posterior` (a numeric vector of log-posterior scores for each iteration). These represent samples drawn from the marginal posterior distribution $P(T,F|G)$ used for Bayesian model averaging.
#'
#' @examples
#' \dontrun{
#' # 1. Define network parameters
#' set.seed(235)
#' num_nodes <- 10
#' sample_size <- 50
#'
#' # 2. Generate true network and simulate data
#' true_network <- GenerateNetwork(num.node = num_nodes)
#'
#' # Set up Beta priors for root-node probabilities and the noise rate
#' prior_para <- matrix(3, nrow = num_nodes + 1, ncol = 2)
#' prior_para[num_nodes + 1, 1] <- 2
#' prior_para[num_nodes + 1, 2] <- 100
#'
#' # Simulate parameters
#' para <- numeric(num_nodes + 1)
#' for (i in 1:(num_nodes + 1)) {
#'   para[i] <- stats::rbeta(1, prior_para[i, 1], prior_para[i, 2])
#' }
#' para[num_nodes + 1] <- 0.1 # Fixed noise rate for simulation
#'
#' error_matrix <- matrix(stats::rbinom(num_nodes * sample_size, 1, para[num_nodes + 1]),
#'   nrow = num_nodes, ncol = sample_size
#' )
#'
#' dummy_data <- GenerateSample(
#'   trans_matrix = true_network,
#'   num.node = num_nodes,
#'   SampleSize = sample_size,
#'   para = para,
#'   error = error_matrix
#' )
#'
#' # 3. Run the MCMC sampler
#' mcmc_results <- run_bbni(
#'   GeneData = dummy_data,
#'   num.node = num_nodes,
#'   SampleSize = sample_size,
#'   prior_para = prior_para,
#'   num_update = 100, # Scaled down for example speed
#'   penalty = 0.1,
#'   prop.ratio = 0.1
#' )
#'
#' # 4. Inspect results
#' tail(mcmc_results$log_posterior)
#' }
#'
#' @importFrom stats runif
#' @export
run_bbni <- function(GeneData, num.node, SampleSize, prior_para,
                     num_update, penalty, prop.ratio) {
  ###############  MCMC
  # prop.ratio<-0.1       # proposal information is used with probability prop.ratio
  prop_beta1 <- 2
  prop_beta2 <- 100
  pseudo.count <- 0.01
  Trans_Func <- seq(1, 12) # all transition functions where each node has two parents at most
  # num.node=20; SampleSize=200 overrides function argument
  # error.prop<-0.2; pseudo.count<-0.01 overrides function argument
  # threshold<-SampleSize*error.prop overrides function argument
  # num_update=200 overrides function argument
  # penalty=0.1 overrides function argument
  # prior_para<-matrix(3, nrow=(num.node+1), ncol=2)
  # prior_para[num.node+1,1]<-2;prior_para[num.node+1,2]<-100    # for error, it should be beta1e<beta2e

  # para<-numeric(); ncp=0
  # for (i in 1:nrow(prior_para))
  # para[i]<-rbeta(1, prior_para[i,1], prior_para[i,2],ncp)
  # para[21]=0.1
  # error<-matrix(0, nrow=num.node,ncol=SampleSize)
  # for(i in 1:num.node)
  #   error[i,]<-rbinom(SampleSize,1,prob=para[num.node+1])
  # ###############
  #   true_network=GenerateNetwork(num.node)       # randomly generate a network

  #   GeneData=matrix(nrow=num.node, ncol=SampleSize)   # based on the generated network, create a data set
  #   GeneData=GenerateSample(true_network)
  # # GeneData<-read.table("GeneData_randomdatanoise0.1proposalinitialpenalty0.1multiplechains_3.txt",header=T)  #true network
  # # GeneData<-data.matrix(GeneData)
  # # true_network<-read.table("TrueNetwork_randomdatanoise0.1proposalinitialpenalty0.1multiplechains_3.txt",header=T)
  # # true_network<-data.matrix(true_network)
  #   true_incid_matrix<-true_network
  #     for (i in 1:nrow(true_incid_matrix))
  #     for (j in 1:ncol(true_incid_matrix))
  #       if (true_incid_matrix[i,j]>0)
  #       true_incid_matrix[i,j]<-1
  #   true_logpost=Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM =true_network)[[1]][length(Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM =true_network))]

  Candidate <- ProposalConstruction(GeneData, SampleSize) # create the proposal for generated data
  prior.triplet <- Candidate[[1]]
  # prior.triplet=read.table("proposal_triplet_true_network.txt",header=T)
  prior.pairwise <- Candidate[[2]]
  # prior.pairwise=read.table("proposal_pairwise_true_network.txt",header=T)
  #####################################
  StandAlone <- matrix(nrow = num.node, ncol = 2) # calculate the BIC for each node
  for (i in 1:num.node)
  {
    succ.prob <- (sum(GeneData[i, ]) + pseudo.count) / (SampleSize + pseudo.count) # calculate BIC for standalone gene
    BIC.value <- -2 * (sum(GeneData[i, ] + pseudo.count) * log(succ.prob) + (SampleSize - sum(GeneData[i, ])) * log(1 - succ.prob)) + 1 * log(SampleSize)
    post.data <- exp(-0.5 * BIC.value)
    StandAlone[i, 1] <- i
    StandAlone[i, 2] <- BIC.value
  }
  # StandAlone=read.table("proposal_standalone_true_network.txt",header=T)
  ###################################  multiple independent chains
  #  All_Trans_Func_Matrix=list()
  #  All_Correct_Rate=list()
  #  All_Logpost=list()
  ###################################

  for (iii in 1:1) # iii: number of simulations
  {
    trans_func_matrix <- ConstructInitial(Candidate, num.node) # use the randomly selected initial
    # trans_func_matrix=read.table("correctrate0.4trans_func_matrix_twostepneighborhood3_8_works_finalnetwork.txt",header=T) # use specific initial
    # trans_func_matrix=read.table("CorrectRate0.55AsStartingNetwork.txt",header=T)
    incid_matrix <- trans_func_matrix
    for (i in 1:nrow(incid_matrix)) {
      for (j in 1:ncol(incid_matrix)) {
        if (incid_matrix[i, j] > 0) {
          incid_matrix[i, j] <- 1
        }
      }
    }
    ances_matrix <- update_ancestor_matrix(incid_matrix)


    Incidence_Matrix <- list()
    Ancestor_Matrix <- list() # Matrix: Ancestor Matrix recording ancesor-offspring relatons for the whole chain
    Trans_Func_Matrix <- list() # Matrix: Transition Function Matrix for the whole chain
    Sample_Matrix <- matrix(nrow = num.node * num_update + 1, ncol = num.node + 2)
    Sample_Matrix[1, ] <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = trans_func_matrix)[[2]]
    Incidence_Matrix[[1]] <- incid_matrix
    Ancestor_Matrix[[1]] <- ances_matrix
    Trans_Func_Matrix[[1]] <- trans_func_matrix
    num <- numeric()
    logpost <- numeric()
    all_logpost <- numeric()
    aa <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = Trans_Func_Matrix[[1]])
    all_logpost[1] <- aa[[1]][length(aa[[1]])]
    logpost[1] <- all_logpost[1]
    n <- 1
    num[1] <- 1
    iter <- 1
    jump_point <- numeric()
    jump_point[1] <- 1 # paramters for each chain
    for (ii in 1:num_update) # run ii full rounds, with each round of num.node times
    {
      if (ii <= round(0.1 * num_update)) { # use adaptive ratio of using proposal information
        prop.ratio <- 0.1
      }
      if (ii > round(0.1 * num_update)) {
        prop.ratio <- 0.9
      }

      update_order <- sample.int(num.node, num.node, replace = FALSE)
      for (k in 1:length(update_order)) #   consider the updating node g_k
      {
        cat("ii", ii, "th iteration is running", "\n")
        cat("n=", n, "\n")
        old <- n
        #   print("iter="); print(iter)
        # cat("true_logpost=", true_logpost, "\n")
        cat("all_logpost=", all_logpost[iter], "\n")
        current_incid_matrix <- Incidence_Matrix[[iter]]
        current_ances_matrix <- Ancestor_Matrix[[iter]]
        current_trans_func_matrix <- Trans_Func_Matrix[[iter]]
        current_post <- all_logpost[iter] # here current_post is a scale

        parent_of_update <- numeric()
        j <- 1 # find the parent for node update_order[k]
        for (i in 1:ncol(current_incid_matrix)) {
          if (current_incid_matrix[update_order[k], i] != 0) {
            parent_of_update[j] <- i
            j <- j + 1
          }
        }
        swap_candi <- numeric()
        j2 <- 1 # swap_candi  is used for swapping parent action     # swap_candi && legal_parent depend on current_parent
        legal_parent <- numeric()
        j1 <- 1 # legal_parent is used for adding parent action
        for (i in 1:num.node) {
          if (i != update_order[k] && current_ances_matrix[i, update_order[k]] != 1 && (i %in% parent_of_update == FALSE)) {
            legal_parent[j1] <- i
            j1 <- j1 + 1 # ; print(i)
          }
        }
        swap_candi <- legal_parent
        #################################    # determine the proposal for current updating node
        pairwise.prior.set <- matrix() # clear variables
        pairwise.prior.pare <- matrix()
        pairwise.prior.set <- prior.pairwise[prior.pairwise[, 2] == update_order[k], ] # pairwise proposal for node update_order[k]
        pairwise.prior.set <- data.matrix(pairwise.prior.set)
        if (ncol(pairwise.prior.set) == 1) {
          pairwise.prior.set <- t(pairwise.prior.set)
          pairwise.prior.pare <- t(data.matrix(pairwise.prior.set[, 1]))
          #     all.maxscore.pairwise=pairwise.prior.set
        }
        if (ncol(pairwise.prior.set) > 1 && nrow(pairwise.prior.set) > 0) {
          pairwise.prior.pare <- pairwise.prior.set[, 1]
          #     all.maxscore.pairwise=pairwise.prior.set[pairwise.prior.set[,4]==max(pairwise.prior.set[,4]),]
        }

        triplet.prior.set <- matrix()
        triplet.prior.pare <- matrix()
        triplet.prior.set <- prior.triplet[prior.triplet[, 3] == update_order[k], ] # triplet proposal for node update_order[k]
        triplet.prior.set <- data.matrix(triplet.prior.set)
        if (ncol(triplet.prior.set) > 1 && nrow(triplet.prior.set) > 0) {
          triplet.prior.pare <- triplet.prior.set[, 1:2]
          #      all.maxscore.triplet=triplet.prior.set[triplet.prior.set[,5]==max(triplet.prior.set[,5]),]
        }
        if (ncol(triplet.prior.set) == 1) {
          triplet.prior.set <- t(triplet.prior.set)
          triplet.prior.pare <- t(data.matrix(triplet.prior.set[, 1:2]))
          #       all.maxscore.triplet=triplet.prior.set
        }
        ###########################  find the minimum mismatches for gene pairs and gene triplet
        #   BIC.pairwise=0; BIC.triplet=0 # case1:  "numeric case"
        #     if(is.numeric(all.maxscore.pairwise)==T && length(all.maxscore.pairwise)>0)
        #        maxscore.pairwise=all.maxscore.pairwise
        #
        #     if(is.numeric(all.maxscore.triplet)==T && length(all.maxscore.triplet)>0)
        #        maxscore.triplet=all.maxscore.triplet
        #
        #     if (is.matrix(all.maxscore.pairwise)==T)
        #      if (nrow(all.maxscore.pairwise)>0)             # randomly select one row if maxscore.pairwsie has multiple rows
        #       if (is.matrix(all.maxscore.triplet)==T)
        #        if(nrow(all.maxscore.triplet)>0)
        #      {
        #       maxscore.pairwise=maxscore.pairwise[sample.int(nrow(maxscore.pairwise),1),]
        #       maxscore.triplet=maxscore.triplet[sample.int(nrow(maxscore.triplet),1),]
        #      }
        #      BIC.pairwise=maxscore.pairwise[5]
        #      BIC.triplet=maxscore.triplet[6]
        #      BIC.standalone=StandAlone[update_order[k],2]
        # pairwise.ratio=BIC.pairwise/(BIC.pairwise+BIC.triplet)
        pairwise.ratio <- 1 / 2
        ############ &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        # case one: no parent
        if (length(parent_of_update) == 0) #  there are only add-parent moves, every time only one move can be proposed
          {
            num_legal_parent <- length(legal_parent)
            if (num_legal_parent > 1) {
              uu <- runif(1)
            }
            if (num_legal_parent == 0) { # if no candiate nodes are available for adding parent, then ignore this node
              uu <- 0
            }
            if (num_legal_parent == 1 || (num_legal_parent > 1 && uu >= pairwise.ratio)) {
              ################################     Every move should have a reasonable proposal probability
              # proposal move 1: add one parent
              prop.legal.overlap <- intersect(pairwise.prior.pare, legal_parent)

              prop.prob <- runif(1)
              aa <- 0
              if (prop.prob >= prop.ratio) { # proposal information is used with probability 1-prop.ratio
                if (length(prop.legal.overlap) >= 1 && nrow(pairwise.prior.set) > 0) {
                  aa <- aa + 1
                  candidate.pare.set <- pairwise.prior.set[pairwise.prior.set[, 1] %in% prop.legal.overlap == T, ]
                  if (is.numeric(candidate.pare.set) == T) {
                    add_parent <- candidate.pare.set
                  }
                  if (is.matrix(candidate.pare.set) == T) { # should modify to use multinomial distribution
                    add_parent <- candidate.pare.set[candidate.pare.set[, 4] == max(candidate.pare.set[, 4]), ]
                  } # use the most likely one
                  if (is.matrix(add_parent) == T) {
                    if (nrow(add_parent) > 1) {
                      add_parent <- add_parent[sample.int(nrow(add_parent), 1), ]
                    }
                  }

                  add_one_parent <- add_parent[1]
                }
              }
              if (prop.prob < prop.ratio || aa == 0) {
                sample_node <- sample.int(length(legal_parent), 1, replace = F)
                add_one_parent <- legal_parent[sample_node]
              }
              prop_incid_matrix <- current_incid_matrix
              prop_incid_matrix[update_order[k], add_one_parent] <- 1
              prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
              prop_trans_func_matrix <- current_trans_func_matrix
              if (aa > 0) {
                prop_trans_func_matrix[update_order[k], add_one_parent] <- add_parent[3]
              } # function prior
              if (aa == 0) {
                func_order <- 10 + sample.int(2, 1)
                prop_trans_func_matrix[update_order[k], add_one_parent] <- func_order
              }
              add_one_prob <- 1 / length(legal_parent)
              prop_trans_func_prob <- 1 / 2 #  there are in all two types of boolean functions to choose for pairwise genes
              xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix)
              prop_sample <- xxx[[2]]
              prop_post <- xxx[[1]][length(xxx[[1]])]
              prop_sample_prob <- numeric()
              p2c_prob <- numeric() # calculate acceptance probability
              curr_sample_prob <- numeric()
              c2p_prob <- numeric()
              prop_sample_prob[1] <- prop_post
              prop_sample_prob[2] <- 0 #  this is the prior && likelihood
              curr_sample_prob[1] <- current_post[1]
              curr_sample_prob[2] <- 0
              prop_sample_prob[3] <- log(1 / 2)
              curr_sample_prob[3] <- log(1) #  this is P(R|T)

              if (num_legal_parent > 1) { # Q(T_c|T_p)
                p2c_prob[1] <- 1 / 3
              } # there are 3 possible moves, add one; remove one; swap one
              if (num_legal_parent == 1) {
                p2c_prob[1] <- 1
              }

              if (num_legal_parent > 1) {
                c2p_prob[1] <- 1 / 2 * add_one_prob
              } # this is Q(T_p|T_c) the probability of move type and specific node involved
              else {
                c2p_prob[1] <- add_one_prob
              }
              p2c_prob[2] <- 1 # Q(R_c|T_c)
              c2p_prob[2] <- prop_trans_func_prob # Q(R_p|T_p)
              p2c_prob[3] <- 1 # Q(\theta_c|T_c, R_c), i.e. the prior probability for new root node
              c2p_prob[3] <- 1 # Q(\theta_p|T_p, R_p)

              nume <- sum((prop_sample_prob)) + sum(log(p2c_prob))
              deno <- sum((curr_sample_prob)) + sum(log(c2p_prob))
              acce_prob <- exp(nume - deno)
              ratio <- runif(1)
              if (ratio <= acce_prob) {
                n <- n + 1
              }
            }
            ##################################
            # proposal move 2 add two parents
            if (num_legal_parent > 1 && uu < pairwise.ratio) # add two parents one time
              {
                prop.prob <- runif(1)
                aa <- 0
                if (prop.prob >= prop.ratio) {
                  if (nrow(triplet.prior.set) >= 1) {
                    candidate.pare.set <- list()
                    for (jj in 1:nrow(triplet.prior.set)) {
                      if (length(intersect(triplet.prior.set[jj, 1:2], legal_parent)) == 2) {
                        aa <- aa + 1
                        candidate.pare.set[[aa]] <- triplet.prior.set[jj, ]
                      }
                    }
                    if (length(candidate.pare.set) > 0) {
                      score <- numeric()
                      for (jj in 1:length(candidate.pare.set)) {
                        score[jj] <- candidate.pare.set[[jj]][5]
                      }
                      bbb <- matrix(nrow = length(candidate.pare.set), ncol = 6)
                      for (jjj in 1:length(candidate.pare.set)) {
                        bbb[jjj, ] <- candidate.pare.set[[jjj]]
                      }

                      max.score.candidate <- bbb[bbb[, 5] == max(score), ]
                      if (is.numeric(max.score.candidate) == T) {
                        add_parent <- max.score.candidate
                      }
                      if (is.matrix(max.score.candidate) == T) {
                        if (nrow(max.score.candidate) > 1) {
                          add_parent <- max.score.candidate[sample.int(nrow(max.score.candidate), 1), ]
                        }
                      }

                      add_two_parent <- c(add_parent[1], add_parent[2])
                    }
                  }
                }
                if (prop.prob < prop.ratio || aa == 0) {
                  sample_two_parent <- sample.int(length(legal_parent), 2, replace = FALSE)
                  add_two_parent <- c(legal_parent[sample_two_parent[1]], legal_parent[sample_two_parent[2]])
                }
                prop_incid_matrix <- current_incid_matrix
                prop_incid_matrix[update_order[k], add_two_parent[1]] <- 1
                prop_incid_matrix[update_order[k], add_two_parent[2]] <- 1
                prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
                prop_trans_func_matrix <- current_trans_func_matrix
                if (aa > 0) {
                  prop_trans_func_matrix[update_order[k], add_two_parent[1]] <- add_parent[4]
                  prop_trans_func_matrix[update_order[k], add_two_parent[2]] <- add_parent[4]
                }
                if (aa == 0) {
                  func_order <- sample.int(10, 1)
                  prop_trans_func_matrix[update_order[k], add_two_parent[1]] <- func_order
                  prop_trans_func_matrix[update_order[k], add_two_parent[2]] <- func_order
                }
                sample_two_prob <- 1 / choose(num_legal_parent, 2)
                prop_trans_func_prob <- 1 / 10
                xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix)
                prop_sample <- xxx[[2]]
                prop_post <- xxx[[1]][length(xxx[[1]])]
                prop_sample_prob <- numeric()
                p2c_prob <- numeric() # calculate acceptance probability
                curr_sample_prob <- numeric()
                c2p_prob <- numeric()
                prop_sample_prob[1] <- prop_post
                prop_sample_prob[2] <- 0 #  this is the prior && likelihood
                curr_sample_prob[1] <- current_post[1]
                curr_sample_prob[2] <- 0
                prop_sample_prob[3] <- log(prop_trans_func_prob)
                curr_sample_prob[3] <- log(1) #  this is P(R|T)

                if (length(swap_candi) > 1) {
                  p2c_prob[1] <- 1 / 4
                } #  this is Q(T_c|T_p) inverse move, because there are 4 moves after adding two parents, depending on the number of nodes to swap
                if (length(swap_candi) == 1) {
                  p2c_prob[1] <- 1 / 3
                }
                if (length(swap_candi) == 0) {
                  p2c_prob[1] <- 1 / 2
                }
                c2p_prob[1] <- 1 / 2 * sample_two_prob # this is Q(T_p|T_c) there are adding one parent, adding two parents, two moves
                p2c_prob[2] <- 1 # Q(R_c|T_c)
                c2p_prob[2] <- prop_trans_func_prob # Q(R_p|T_p)
                p2c_prob[3] <- 1 # Q(\theta_c|T_c, R_c), i.e. the prior probability for new root node
                c2p_prob[3] <- 1 # Q(\theta_p|T_p, R_p)

                nume <- sum((prop_sample_prob)) + sum(log(p2c_prob))
                deno <- sum((curr_sample_prob)) + sum(log(c2p_prob))
                acce_prob <- exp(nume - deno)
                ratio <- runif(1)
                if (ratio <= acce_prob) {
                  n <- n + 1
                }
              }
          }
        ########## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        # case 2: one parent
        if (length(parent_of_update) == 1) {
          uu <- runif(1)
          # one.parent.ratio=2/5*(1/BIC.triplet/(1/BIC.triplet+1/BIC.standalone))
          one.parent.ratio <- 1 / 6
          #################################
          # proposal move 1: add one parent
          if (length(legal_parent) > 0 && uu < one.parent.ratio) {
            prop.prob <- runif(1)
            aa <- 0
            if (prop.prob >= prop.ratio) {
              if (nrow(triplet.prior.set) >= 1 && ncol(triplet.prior.set) > 1) {
                candidate.pare.set <- list()
                for (jj in 1:nrow(triplet.prior.set)) {
                  if (length(intersect(parent_of_update, triplet.prior.pare[jj, ])) == 1) {
                    if (parent_of_update %in% triplet.prior.pare[jj, ] == T && setdiff(triplet.prior.pare[jj, ], parent_of_update) %in% legal_parent == T) {
                      aa <- aa + 1
                      candidate.pare.set[[aa]] <- triplet.prior.set[jj, ]
                    }
                  }
                }
                if (length(candidate.pare.set) > 0) {
                  score <- numeric()
                  for (jj in 1:length(candidate.pare.set)) {
                    score[jj] <- candidate.pare.set[[jj]][5]
                  }
                  bbb <- matrix(nrow = length(candidate.pare.set), ncol = 6)
                  for (jjj in 1:length(candidate.pare.set)) {
                    bbb[jjj, ] <- candidate.pare.set[[jjj]]
                  }

                  max.score.candidate <- bbb[bbb[, 5] == max(score), ]
                  if (is.numeric(max.score.candidate) == T) {
                    add_parent <- max.score.candidate
                  }
                  if (is.matrix(max.score.candidate) == T) {
                    if (nrow(max.score.candidate) > 1) {
                      add_parent <- max.score.candidate[sample.int(nrow(max.score.candidate), 1), ]
                    }
                  }

                  add_one_parent <- setdiff(add_parent[1:2], parent_of_update)
                }
              }
            }
            if (prop.prob < prop.ratio || aa == 0) {
              sample_node <- sample.int(length(legal_parent), 1, replace = FALSE)
              add_one_parent <- legal_parent[sample_node]
            }
            prop_incid_matrix <- current_incid_matrix
            prop_incid_matrix[update_order[k], add_one_parent] <- 1
            prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
            prop_trans_func_matrix <- current_trans_func_matrix
            if (aa > 0) {
              prop_trans_func_matrix[update_order[k], add_one_parent] <- add_parent[4]
              prop_trans_func_matrix[update_order[k], parent_of_update] <- add_parent[4]
            }
            if (aa == 0) {
              func_order <- sample.int(10, 1)
              prop_trans_func_matrix[update_order[k], add_one_parent] <- func_order
              prop_trans_func_matrix[update_order[k], parent_of_update] <- func_order
            }
            add_one_prob <- 1 / length(legal_parent)
            prop_trans_func_prob <- 1 / 10 # there are 10 functions to choose
            xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix)
            prop_sample <- xxx[[2]]
            prop_post <- xxx[[1]][length(xxx[[1]])]
            prop_sample_prob <- numeric()
            p2c_prob <- numeric() # calculate acceptance probability
            curr_sample_prob <- numeric()
            c2p_prob <- numeric()
            prop_sample_prob[1] <- prop_post
            prop_sample_prob[2] <- 0 #  this is the prior && likelihood
            curr_sample_prob[1] <- current_post[1]
            curr_sample_prob[2] <- 0
            prop_sample_prob[3] <- log(prop_trans_func_prob)
            curr_sample_prob[3] <- log(1) #  this is P(R|T)

            if (length(swap_candi) > 1) {
              p2c_prob[1] <- 1 / 4
            } #  this is Q(T_c|T_p) inverse move, because there are 4 moves after adding one parents, depending on the number of nodes to swap
            if (length(swap_candi) == 1) {
              p2c_prob[1] <- 1 / 3
            }
            if (length(swap_candi) == 0) {
              p2c_prob[1] <- 1 / 2
            }

            if (length(swap_candi) > 0) {
              c2p_prob[1] <- 1 / 3 * add_one_prob
            } # this is Q(T_p|T_c) probability of which move and which node selected
            if (length(swap_candi) == 0) {
              c2p_prob[1] <- 1 / 2 * add_one_prob
            }

            p2c_prob[2] <- 1 # Q(R_c|T_c)
            c2p_prob[2] <- prop_trans_func_prob # Q(R_p|T_p)
            p2c_prob[3] <- 1 # Q(\theta_c|T_c, R_c), i.e. the prior probability for new root node
            c2p_prob[3] <- 1 # Q(\theta_p|T_p, R_p)

            nume <- sum((prop_sample_prob)) + sum(log(p2c_prob))
            deno <- sum((curr_sample_prob)) + sum(log(c2p_prob))
            acce_prob <- exp(nume - deno)
            ratio <- runif(1)
            if (ratio <= acce_prob) {
              n <- n + 1
            }
          }
          ##################################
          # proposal move 2, swap one parent
          if (length(swap_candi) > 0 && uu >= 1 / 6 && uu < 2 / 6) {
            swap.prior.overlap <- intersect(swap_candi, pairwise.prior.pare)
            prop.prob <- runif(1)
            aa <- 0
            if (prop.prob > prop.ratio) {
              if (length(swap.prior.overlap) > 0 && nrow(pairwise.prior.set) > 0) {
                aa <- aa + 1
                swap.candidate <- pairwise.prior.set[pairwise.prior.set[, 1] %in% swap.prior.overlap == T, ]
                if (is.numeric(swap.candidate) == T) { # i.e. swap.candidate has one row
                  swap_parent <- swap.candidate
                }
                if (is.matrix(swap.candidate) == T) {
                  swap_parent <- swap.candidate[swap.candidate[, 4] == max(swap.candidate[, 4]), ]
                } # multiple norw may have same weight

                if (is.matrix(swap_parent) == T) {
                  if (nrow(swap_parent) > 1) {
                    swap_parent <- swap_parent[sample.int(nrow(swap_parent), 1), ]
                  }
                } # every candidate has equal chance
                swap_one_node <- swap_parent[1]
              }
            }
            if (prop.prob < prop.ratio || aa == 0) {
              sample_node <- sample.int(length(swap_candi), 1, replace = FALSE)
              swap_one_node <- swap_candi[sample_node]
            }
            prop_incid_matrix <- current_incid_matrix
            prop_incid_matrix[update_order[k], parent_of_update] <- 0
            prop_incid_matrix[update_order[k], swap_one_node] <- 1
            prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
            prop_trans_func_matrix <- current_trans_func_matrix
            prop_trans_func_matrix[update_order[k], parent_of_update] <- 0
            if (aa > 0) {
              prop_trans_func_matrix[update_order[k], swap_one_node] <- swap_parent[3]
            }
            if (aa == 0) {
              func_order <- 10 + sample.int(2, 1)
              prop_trans_func_matrix[update_order[k], swap_one_node] <- func_order
            }
            swap_one_prob <- 1 / length(swap_candi)
            prop_trans_func_prob <- 1 / 2 #  new proposal transition function matrix
            xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix)
            prop_sample <- xxx[[2]]
            prop_post <- xxx[[1]][length(xxx[[1]])]
            prop_sample_prob <- numeric()
            p2c_prob <- numeric() # calculate acceptance probability
            curr_sample_prob <- numeric()
            c2p_prob <- numeric()
            prop_sample_prob[1] <- prop_post
            prop_sample_prob[2] <- 0 #  this is the prior && likelihood
            curr_sample_prob[1] <- current_post[1]
            curr_sample_prob[2] <- 0
            prop_sample_prob[3] <- log(prop_trans_func_prob)
            curr_sample_prob[3] <- log(1) #  this is P(R|T)

            if (length(legal_parent) > 0) {
              p2c_prob[1] <- 1 / 3
            } #  this is Q(T_c|T_p) inverse move, because there are 4 moves after adding one parents, depending on the number of nodes to swap
            if (length(legal_parent) == 0) {
              p2c_prob[1] <- 1 / 2
            }

            if (length(legal_parent) > 0) {
              c2p_prob[1] <- 1 / 3 * swap_one_prob
            } # this is Q(T_p|T_c) probability of which move and which node selected
            if (length(legal_parent) == 0) {
              c2p_prob[1] <- 1 / 2 * swap_one_prob
            }
            p2c_prob[2] <- 1 # Q(R_c|T_c)
            c2p_prob[2] <- prop_trans_func_prob # Q(R_p|T_p)
            p2c_prob[3] <- 1 # Q(\theta_c|T_c, R_c) no root node is introduced.
            c2p_prob[3] <- 1 # Q(\theta_p|T_p, R_p)

            nume <- sum((prop_sample_prob)) + sum(log(p2c_prob))
            deno <- sum((curr_sample_prob)) + sum(log(c2p_prob))
            acce_prob <- exp(nume - deno)
            ratio <- runif(1)
            if (ratio <= acce_prob) {
              n <- n + 1
            }
          }
          ##################################
          # proposal move 3; remove one parent ### this move may introduce new root node
          if (uu >= 2 / 6 && uu < 3 / 6) {
            remove_one_node <- parent_of_update
            remove_one_prob <- 1
            prop_incid_matrix <- current_incid_matrix
            prop_incid_matrix[update_order[k], remove_one_node] <- 0 # remove parent move does not need function.
            prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
            prop_trans_func_matrix <- current_trans_func_matrix
            prop_trans_func_matrix[update_order[k], remove_one_node] <- 0
            prop_trans_func_prob <- 1
            xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix)
            prop_sample <- xxx[[2]]
            prop_post <- xxx[[1]][length(xxx[[1]])]
            prop_sample_prob <- numeric()
            p2c_prob <- numeric() # calculate acceptance probability
            curr_sample_prob <- numeric()
            c2p_prob <- numeric()
            prop_sample_prob[1] <- prop_post
            prop_sample_prob[2] <- 0 #  this is the prior && likelihood
            curr_sample_prob[1] <- current_post[1]
            curr_sample_prob[2] <- 0
            prop_sample_prob[3] <- log(prop_trans_func_prob)
            curr_sample_prob[3] <- log(1) #  this is P(R|T)

            if (length(legal_parent) == 0) {
              p2c_prob[1] <- 1
            }
            if (length(legal_parent) == 1) {
              p2c_prob[1] <- 1
            } #  this is Q(T_c|T_p) inverse move
            if (length(legal_parent) > 1) {
              p2c_prob[1] <- 1 / 2
            }

            if (length(legal_parent) > 0) {
              c2p_prob[1] <- 1 / 3 * remove_one_prob
            } # this is Q(T_p|T_c) probability of which move and which node selected
            if (length(legal_parent) == 0) {
              c2p_prob[1] <- 1 / 2 * remove_one_prob
            }
            p2c_prob[2] <- 1 / 2 # Q(R_c|T_c)
            c2p_prob[2] <- prop_trans_func_prob # Q(R_p|T_p)
            p2c_prob[3] <- 1 # Q(\theta_c|T_c, R_c)
            c2p_prob[3] <- 1 # Q(\theta_p|T_p, R_p)

            nume <- sum((prop_sample_prob)) + sum(log(p2c_prob))
            deno <- sum((curr_sample_prob)) + sum(log(c2p_prob))
            acce_prob <- exp(nume - deno)
            ratio <- runif(1)
            if (ratio <= acce_prob) {
              n <- n + 1
            }
          }
          ####################################
          # additional move 1: reverse one arc in pairwise genes
          if (uu >= 3 / 6 && uu < 4 / 6) {
            prop_incid_matrix <- current_incid_matrix
            prop_incid_matrix[update_order[k], parent_of_update] <- 0
            prop_incid_matrix[parent_of_update, update_order[k]] <- 1
            parent_parent <- numeric()
            j <- 0
            for (i in 1:nrow(current_incid_matrix)) {
              if (current_incid_matrix[parent_of_update, i] > 0) {
                j <- j + 1
                parent_parent[j] <- i
              }
            }
            prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
            if (check_ances_matrix(prop_ances_matrix) == 0 && length(parent_parent) == 0) # make sure no loops && parent_of_update no parents
              {
                prop_trans_func_matrix <- current_trans_func_matrix
                func_order <- current_trans_func_matrix[update_order[k], parent_of_update]
                prop_trans_func_matrix[update_order[k], parent_of_update] <- 0
                prop_trans_func_matrix[parent_of_update, update_order[k]] <- func_order
                xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix)
                prop_sample <- xxx[[2]]
                prop_post <- xxx[[1]][length(xxx[[1]])]

                nume <- sum((prop_post))
                deno <- sum((current_post))
                acce_prob <- exp(nume - deno)
                ratio <- runif(1)
                if (ratio <= acce_prob) {
                  n <- n + 1
                }
              }
          }
          ###################################
          # additional move 2: reverse two arcs simultanously
          if (uu >= 4 / 6 && uu < 5 / 6) {
            parent_parent <- numeric()
            j1 <- 0
            children <- numeric()
            j <- 0 # find the children of parent_of_update
            for (i in 1:ncol(current_incid_matrix))
            {
              if (i != update_order[k] && current_incid_matrix[i, parent_of_update] > 0) {
                j <- j + 1
                children[j] <- i
              }
              if (current_incid_matrix[parent_of_update, i] > 0) {
                j1 <- j1 + 1
                parent_parent[j1] <- i
              }
            }
            if (length(children) == 1 && length(parent_parent) == 0) # parent_of_update has two exact 2 children and no parents
              {
                parent_child <- numeric()
                j <- 0 # find parent of children to determine the trans_func_matrix
                for (i in 1:num.node) {
                  if (current_incid_matrix[children, i] > 0) {
                    j <- j + 1
                    parent_child[j] <- i
                  }
                }

                prop_incid_matrix <- current_incid_matrix
                prop_incid_matrix[update_order[k], parent_of_update] <- 0
                prop_incid_matrix[children, parent_of_update] <- 0
                prop_incid_matrix[parent_of_update, update_order[k]] <- 1
                prop_incid_matrix[parent_of_update, children] <- 1
                prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
                if (check_ances_matrix(prop_ances_matrix) == 0) # ensure no directed cycles
                  {
                    prop_trans_func_matrix <- current_trans_func_matrix
                    prop_trans_func_matrix[update_order[k], parent_of_update] <- 0
                    prop_trans_func_matrix[children, parent_of_update] <- 0
                    func_order <- sample.int(10, 1)
                    prop_trans_func_matrix[parent_of_update, update_order[k]] <- func_order
                    prop_trans_func_matrix[parent_of_update, children] <- func_order
                    if (length(parent_child) == 2) {
                      prop_trans_func_matrix[children, setdiff(parent_child, parent_of_update)] <- 10 + sample.int(2, 1)
                    }
                    xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix)
                    prop_sample <- xxx[[2]]
                    prop_post <- xxx[[1]][length(xxx[[1]])]

                    nume <- sum((prop_post))
                    deno <- sum((current_post))
                    acce_prob <- exp(nume - deno)
                    ratio <- runif(1)
                    if (ratio <= acce_prob) {
                      n <- n + 1
                    }
                  }
              }
          }
          ############################################# additional move: swap two with one
          if (uu >= 5 / 6 && length(swap_candi) > 1) {
            legal.triplet.pare <- intersect(swap_candi, triplet.prior.pare)
            prop.prob <- runif(1)
            aa <- 0
            if (prop.prob > prop.ratio) {
              swap.candidate <- list()
              if (length(intersect(parent_of_update, legal.triplet.pare)) == 0 && nrow(triplet.prior.set) > 0) {
                for (jj in 1:nrow(triplet.prior.pare)) {
                  if (length(intersect(parent_of_update, triplet.prior.pare[jj, ])) == 0 && length(intersect(triplet.prior.pare[jj, ], swap_candi)) == 2) {
                    aa <- aa + 1
                    swap.candidate[[aa]] <- triplet.prior.set[jj, ]
                  }
                }
              }
              if (length(swap.candidate) > 0) {
                score <- numeric()
                for (jj in 1:length(swap.candidate)) {
                  score[jj] <- swap.candidate[[jj]][5]
                }
                bbb <- matrix(nrow = length(swap.candidate), ncol = 6) # convert list to matrix
                for (jjj in 1:length(swap.candidate)) {
                  bbb[jjj, ] <- swap.candidate[[jjj]]
                }
                max.score.candidate <- bbb[bbb[, 5] == max(score), ]
                if (is.numeric(max.score.candidate) == T) {
                  swap_parent <- max.score.candidate
                }
                if (is.matrix(max.score.candidate) == T) {
                  if (nrow(max.score.candidate) > 1) {
                    swap_parent <- max.score.candidate[sample.int(nrow(max.score.candidate), 1), ]
                  }
                }
                swap_two_parent <- swap_parent[1:2]
              }
            }
            if (prop.prob < prop.ratio || aa == 0) {
              sample_two_node <- sample.int(length(swap_candi), 2, replace = FALSE)
              swap_two_parent <- c(swap_candi[sample_two_node[1]], swap_candi[sample_two_node[2]])
            }
            prop_incid_matrix <- current_incid_matrix
            prop_incid_matrix[update_order[k], parent_of_update] <- 0
            prop_incid_matrix[update_order[k], swap_two_parent[1]] <- 1
            prop_incid_matrix[update_order[k], swap_two_parent[2]] <- 1
            prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
            prop_trans_func_matrix <- current_trans_func_matrix
            prop_trans_func_matrix[update_order[k], parent_of_update] <- 0
            if (aa > 0) {
              prop_trans_func_matrix[update_order[k], swap_two_parent[1]] <- swap_parent[4]
              prop_trans_func_matrix[update_order[k], swap_two_parent[2]] <- swap_parent[4]
            }
            if (aa == 0) {
              func_order <- sample.int(10, 1)
              prop_trans_func_matrix[update_order[k], swap_two_parent[1]] <- func_order
              prop_trans_func_matrix[update_order[k], swap_two_parent[2]] <- func_order
            }
            swap_two_prob <- 1 / choose(length(swap_candi), 2)
            prop_trans_func_prob <- 1 / 10
            xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix)
            prop_sample <- xxx[[2]]
            prop_post <- xxx[[1]][length(xxx[[1]])]
            prop_sample_prob <- numeric()
            p2c_prob <- numeric() # calculate acceptance probability
            curr_sample_prob <- numeric()
            c2p_prob <- numeric()

            prop_sample_prob[1] <- prop_post
            prop_sample_prob[2] <- 0 #  this is the prior && likelihood
            curr_sample_prob[1] <- current_post[1]
            curr_sample_prob[2] <- 0
            prop_sample_prob[3] <- log(prop_trans_func_prob)
            curr_sample_prob[3] <- log(1) #  this is P(R|T)

            if (length(swap_candi) > 1) {
              p2c_prob[1] <- 1 / 4
            } #  this is Q(T_c|T_p) inverse move, because there are 4 moves after adding one parents, depending on the number of nodes to swap
            if (length(swap_candi) == 1) {
              p2c_prob[1] <- 1 / 3
            }

            if (length(swap_candi) > 1) {
              c2p_prob[1] <- 1 / 4 * swap_two_prob
            } # this is Q(T_p|T_c) probability of which move and which node selected
            if (length(swap_candi) == 1) {
              c2p_prob[1] <- 1 / 3 * swap_two_prob
            }
            p2c_prob[2] <- 1 / 10 # Q(R_c|T_c)
            c2p_prob[2] <- prop_trans_func_prob # Q(R_p|T_p)
            p2c_prob[3] <- 1 # Q(\theta_c|T_c, R_c) no root node is introduced.
            c2p_prob[3] <- 1 # Q(\theta_p|T_p, R_p)

            nume <- sum((prop_sample_prob)) + sum(log(p2c_prob))
            deno <- sum((curr_sample_prob)) + sum(log(c2p_prob))
            acce_prob <- exp(nume - deno)
            ratio <- runif(1)
            if (ratio <= acce_prob) {
              n <- n + 1
            }
          }
        }
        ######### &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        # case3: two parents
        if (length(parent_of_update) == 2) {
          uu <- runif(1) # two.parent.ratio=2/7*(1/BIC.pairwise/(1/BIC.pairwise+1/BIC.standalone))
          two.parent.ratio <- 1 / 7
          ##########################
          # proposal move 1: remove one parent
          if (uu < two.parent.ratio) {
            pare.prior.overlap <- intersect(parent_of_update, pairwise.prior.pare)
            prop.prob <- runif(1)
            aa <- 0
            if (prop.prob >= prop.ratio) {
              if (length(pare.prior.overlap) == 1 && nrow(pairwise.prior.set) > 0) {
                aa <- aa + 1
                remove_one_node <- setdiff(parent_of_update, pairwise.prior.pare)
                remove_pare <- pairwise.prior.set[pairwise.prior.set[, 1] == setdiff(parent_of_update, remove_one_node), ]
                if (is.numeric(remove_pare) == T) {
                  func_order <- remove_pare[3]
                }
                if (is.matrix(remove_pare) == T) {
                  Remove_Pare <- remove_pare[remove_pare[, 4] == max(remove_pare[, 4]), ]
                  if (is.matrix(Remove_Pare) == T) {
                    if (nrow(Remove_Pare) > 1) {
                      Remove_Pare <- Remove_Pare[sample.int(nrow(Remove_Pare), 1), ]
                    }
                  } # each candidate has equal chance
                  func_order <- Remove_Pare[3]
                }
              }
            }
            if (prop.prob < prop.ratio || aa == 0) {
              sample_node <- sample.int(length(parent_of_update), 1, replace = FALSE)
              remove_one_node <- parent_of_update[sample_node]
            }
            remove_one_prob <- 1 / length(parent_of_update)
            prop_incid_matrix <- current_incid_matrix
            prop_incid_matrix[update_order[k], remove_one_node] <- 0
            prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
            prop_trans_func_matrix <- current_trans_func_matrix
            prop_trans_func_matrix[update_order[k], remove_one_node] <- 0
            if (aa > 0) {
              prop_trans_func_matrix[update_order[k], setdiff(parent_of_update, remove_one_node)] <- func_order
            } # by the definition
            if (aa == 0) {
              func_order <- 10 + sample.int(2, 1)
            }
            prop_trans_func_matrix[update_order[k], setdiff(parent_of_update, remove_one_node)] <- func_order
            prop_trans_func_prob <- 1 / 2
            xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix)
            prop_sample <- xxx[[2]]
            prop_post <- xxx[[1]][length(xxx[[1]])]
            prop_sample_prob <- numeric()
            p2c_prob <- numeric() # calculate acceptance probability
            curr_sample_prob <- numeric()
            c2p_prob <- numeric()

            prop_sample_prob[1] <- prop_post
            prop_sample_prob[2] <- 0 #  this is the prior && likelihood
            curr_sample_prob[1] <- current_post[1]
            curr_sample_prob[2] <- 0
            prop_sample_prob[3] <- log(prop_trans_func_prob)
            curr_sample_prob[3] <- log(1) #  this is P(R|T)

            if (length(swap_candi) > 0) {
              p2c_prob[1] <- 1 / 3
            } #  this is Q(T_c|T_p) inverse move, because there are 4 moves after adding one parents, depending on the number of nodes to swap
            if (length(swap_candi) == 0) {
              p2c_prob[1] <- 1 / 2
            }

            if (length(swap_candi) == 0) {
              c2p_prob[1] <- 1 / 2
            }
            if (length(swap_candi) > 1) {
              c2p_prob[1] <- 1 / 4 * remove_one_prob
            } # this is Q(T_p|T_c) probability of which move and which node selected
            if (length(swap_candi) == 1) {
              c2p_prob[1] <- 1 / 3 * remove_one_prob
            }

            p2c_prob[2] <- 1 / 10 # Q(R_c|T_c)
            c2p_prob[2] <- prop_trans_func_prob # Q(R_p|T_p)
            p2c_prob[3] <- 1 # Q(\theta_c|T_c, R_c) no root node is introduced.
            c2p_prob[3] <- 1 # Q(\theta_p|T_p, R_p)

            nume <- sum((prop_sample_prob)) + sum(log(p2c_prob))
            deno <- sum((curr_sample_prob)) + sum(log(c2p_prob))
            acce_prob <- exp(nume - deno)
            ratio <- runif(1)
            if (ratio <= acce_prob) {
              n <- n + 1
            }
          }
          ############################
          # proposal move 2; remove two parents
          if (uu >= two.parent.ratio && uu < 2 / 7) {
            remove_node <- parent_of_update
            remove_prob <- 1
            prop_incid_matrix <- current_incid_matrix
            prop_incid_matrix[update_order[k], remove_node[1]] <- 0
            prop_incid_matrix[update_order[k], remove_node[2]] <- 0
            prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
            prop_trans_func_matrix <- current_trans_func_matrix
            prop_trans_func_matrix[update_order[k], remove_node[1]] <- 0
            prop_trans_func_matrix[update_order[k], remove_node[2]] <- 0
            prop_trans_func_prob <- 1
            xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix)
            prop_sample <- xxx[[2]]
            prop_post <- xxx[[1]][length(xxx[[1]])]
            prop_sample_prob <- numeric()
            p2c_prob <- numeric() # calculate acceptance probability
            curr_sample_prob <- numeric()
            c2p_prob <- numeric()
            prop_sample_prob[1] <- prop_post
            prop_sample_prob[2] <- 0 #  this is the prior && likelihood
            curr_sample_prob[1] <- current_post[1]
            curr_sample_prob[2] <- 0
            prop_sample_prob[3] <- log(prop_trans_func_prob)
            curr_sample_prob[3] <- log(1) #  this is P(R|T)

            p2c_prob[1] <- 1 / 2 #  this is Q(T_c|T_p) inverse move
            c2p_prob[1] <- 1 / 4 * remove_prob # this is Q(T_p|T_c) probability of which move and which node selected
            p2c_prob[2] <- 1 / 10 # Q(R_c|T_c)
            c2p_prob[2] <- prop_trans_func_prob # Q(R_p|T_p)
            p2c_prob[3] <- 1 # Q(\theta_c|T_c, R_c) no root node is introduced.
            c2p_prob[3] <- 1 # Q(\theta_p|T_p, R_p)

            nume <- sum((prop_sample_prob)) + sum(log(p2c_prob))
            deno <- sum((curr_sample_prob)) + sum(log(c2p_prob))
            acce_prob <- exp(nume - deno)
            ratio <- runif(1)
            if (ratio <= acce_prob) {
              n <- n + 1
            }
          }
          ############################
          # proposal move 3: swap one parent
          if (uu >= 2 / 7 && uu < 3 / 7 && length(swap_candi) > 0) {
            legal.triplet.pare <- intersect(swap_candi, triplet.prior.pare)
            prop.prob <- runif(1)
            aa <- 0
            if (prop.prob >= prop.ratio) {
              swap.candidate <- list()
              if (length(intersect(parent_of_update, legal.triplet.pare)) == 1 && nrow(triplet.prior.set) > 0) {
                for (jj in 1:nrow(triplet.prior.pare)) {
                  if (length(intersect(parent_of_update, triplet.prior.pare[jj, ])) == 1) {
                    aa <- aa + 1
                    swap.candidate[[aa]] <- triplet.prior.set[jj, ]
                  }
                }
              }
              if (length(swap.candidate) > 0) {
                score <- numeric()
                for (jj in 1:length(swap.candidate)) {
                  score[jj] <- swap.candidate[[jj]][5]
                }
                bbb <- matrix(nrow = length(swap.candidate), ncol = 6)
                for (jjj in 1:length(swap.candidate)) {
                  bbb[jjj, ] <- swap.candidate[[jjj]]
                }

                max.score.candidate <- bbb[bbb[, 5] == max(score), ]
                if (is.numeric(max.score.candidate) == T) {
                  swap_parent <- max.score.candidate
                }
                if (is.matrix(max.score.candidate) == T) {
                  if (nrow(max.score.candidate) > 1) {
                    swap_parent <- max.score.candidate[sample.int(nrow(max.score.candidate), 1), ]
                  }
                }
                swap_one_node <- setdiff(swap_parent[1:2], parent_of_update)
                sample_one_parent <- setdiff(parent_of_update, swap_parent[1:2])
              }
            }
            if (prop.prob < prop.ratio || aa == 0) {
              sample_node <- sample.int(length(swap_candi), 1, replace = FALSE)
              sample_parent <- sample.int(2, 1, replace = FALSE)
              swap_one_node <- swap_candi[sample_node]
              sample_one_parent <- parent_of_update[sample_parent]
            }
            prop_incid_matrix <- current_incid_matrix
            prop_incid_matrix[update_order[k], sample_one_parent] <- 0
            prop_incid_matrix[update_order[k], swap_one_node] <- 1
            prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
            prop_trans_func_matrix <- current_trans_func_matrix
            prop_trans_func_matrix[update_order[k], sample_one_parent] <- 0
            if (aa > 0) {
              prop_trans_func_matrix[update_order[k], setdiff(parent_of_update, sample_one_parent)] <- swap_parent[4]
              prop_trans_func_matrix[update_order[k], swap_one_node] <- swap_parent[4]
            }
            if (aa == 0) {
              func_order <- sample.int(10, 1)
              prop_trans_func_matrix[update_order[k], setdiff(parent_of_update, sample_one_parent)] <- func_order
              prop_trans_func_matrix[update_order[k], swap_one_node] <- func_order
            }
            swap_one_prob <- 1 / length(swap_candi) * 1 / 2
            prop_trans_func_prob <- 1 / 10
            xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix)
            prop_sample <- xxx[[2]]
            prop_post <- xxx[[1]][length(xxx[[1]])]
            prop_sample_prob <- numeric()
            p2c_prob <- numeric() # calculate acceptance probability
            curr_sample_prob <- numeric()
            c2p_prob <- numeric()
            prop_sample_prob[1] <- prop_post
            prop_sample_prob[2] <- 0 #  this is the prior && likelihood
            curr_sample_prob[1] <- current_post[1]
            curr_sample_prob[2] <- 0
            prop_sample_prob[3] <- log(prop_trans_func_prob)
            curr_sample_prob[3] <- log(1) #  this is P(R|T)

            if (length(swap_candi) == 0) {
              p2c_prob[1] <- 1
            } # since we use log(p2c_prob)
            if (length(swap_candi) > 1) {
              p2c_prob[1] <- 1 / 4
            } #  this is Q(T_c|T_p) inverse move, because there are 4 moves after adding one parents, depending on the number of nodes to swap
            if (length(swap_candi) == 1) {
              p2c_prob[1] <- 1 / 3
            }

            if (length(swap_candi) == 0) {
              c2p_prob[1] <- 1
            }
            if (length(swap_candi) > 1) {
              c2p_prob[1] <- 1 / 4 * swap_one_prob
            } # this is Q(T_p|T_c) probability of which move and which node selected
            if (length(swap_candi) == 1) {
              c2p_prob[1] <- 1 / 3 * swap_one_prob
            }
            p2c_prob[2] <- 1 / 10 # Q(R_c|T_c)
            c2p_prob[2] <- prop_trans_func_prob # Q(R_p|T_p)
            p2c_prob[3] <- 1 # Q(\theta_c|T_c, R_c) no root node is introduced.
            c2p_prob[3] <- 1 # Q(\theta_p|T_p, R_p)

            nume <- sum((prop_sample_prob)) + sum(log(p2c_prob))
            deno <- sum((curr_sample_prob)) + sum(log(c2p_prob))
            acce_prob <- exp(nume - deno)
            ratio <- runif(1)
            if (ratio <= acce_prob) {
              n <- n + 1
            }
          }
          ##############################
          # proposal move 4, swap two parents
          if (uu >= 3 / 7 && uu < 4 / 7 && length(swap_candi) > 1) {
            legal.triplet.pare <- intersect(swap_candi, triplet.prior.pare)
            prop.prob <- runif(1)
            aa <- 0
            if (prop.prob > prop.ratio) {
              swap.candidate <- list()
              if (length(intersect(parent_of_update, legal.triplet.pare)) == 0 && nrow(triplet.prior.set) > 0) {
                for (jj in 1:nrow(triplet.prior.pare)) {
                  if (length(intersect(parent_of_update, triplet.prior.pare[jj, ])) == 0 && length(intersect(triplet.prior.pare[jj, ], swap_candi)) == 2) {
                    aa <- aa + 1
                    swap.candidate[[aa]] <- triplet.prior.set[jj, ]
                  }
                }
              }
              if (length(swap.candidate) > 0) {
                score <- numeric()
                for (jj in 1:length(swap.candidate)) {
                  score[jj] <- swap.candidate[[jj]][5]
                }
                bbb <- matrix(nrow = length(swap.candidate), ncol = 6) # convert list to matrix
                for (jjj in 1:length(swap.candidate)) {
                  bbb[jjj, ] <- swap.candidate[[jjj]]
                }
                max.score.candidate <- bbb[bbb[, 5] == max(score), ]
                if (is.numeric(max.score.candidate) == T) {
                  swap_parent <- max.score.candidate
                }
                if (is.matrix(max.score.candidate) == T) {
                  if (nrow(max.score.candidate) > 1) {
                    swap_parent <- max.score.candidate[sample.int(nrow(max.score.candidate), 1), ]
                  }
                }
                swap_two_parent <- swap_parent[1:2]
              }
            }
            if (prop.prob < prop.ratio || aa == 0) {
              sample_two_node <- sample.int(length(swap_candi), 2, replace = FALSE)
              swap_two_parent <- c(swap_candi[sample_two_node[1]], swap_candi[sample_two_node[2]])
            }
            prop_incid_matrix <- current_incid_matrix
            prop_incid_matrix[update_order[k], parent_of_update[1]] <- 0
            prop_incid_matrix[update_order[k], parent_of_update[2]] <- 0
            prop_incid_matrix[update_order[k], swap_two_parent[1]] <- 1
            prop_incid_matrix[update_order[k], swap_two_parent[2]] <- 1
            prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
            prop_trans_func_matrix <- current_trans_func_matrix
            prop_trans_func_matrix[update_order[k], parent_of_update[1]] <- 0
            prop_trans_func_matrix[update_order[k], parent_of_update[2]] <- 0
            if (aa > 0) {
              prop_trans_func_matrix[update_order[k], swap_two_parent[1]] <- swap_parent[4]
              prop_trans_func_matrix[update_order[k], swap_two_parent[2]] <- swap_parent[4]
            }
            if (aa == 0) {
              func_order <- sample.int(10, 1)
              prop_trans_func_matrix[update_order[k], swap_two_parent[1]] <- func_order
              prop_trans_func_matrix[update_order[k], swap_two_parent[2]] <- func_order
            }
            swap_two_prob <- 1 / choose(length(swap_candi), 2)
            prop_trans_func_prob <- 1 / 10
            xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix)
            prop_sample <- xxx[[2]]
            prop_post <- xxx[[1]][length(xxx[[1]])]
            prop_sample_prob <- numeric()
            p2c_prob <- numeric() # calculate acceptance probability
            curr_sample_prob <- numeric()
            c2p_prob <- numeric()

            prop_sample_prob[1] <- prop_post
            prop_sample_prob[2] <- 0 #  this is the prior && likelihood
            curr_sample_prob[1] <- current_post[1]
            curr_sample_prob[2] <- 0
            prop_sample_prob[3] <- log(prop_trans_func_prob)
            curr_sample_prob[3] <- log(1) #  this is P(R|T)

            if (length(swap_candi) > 1) {
              p2c_prob[1] <- 1 / 4
            } #  this is Q(T_c|T_p) inverse move, because there are 4 moves after adding one parents, depending on the number of nodes to swap
            if (length(swap_candi) == 1) {
              p2c_prob[1] <- 1 / 3
            }

            if (length(swap_candi) > 1) {
              c2p_prob[1] <- 1 / 4 * swap_two_prob
            } # this is Q(T_p|T_c) probability of which move and which node selected
            if (length(swap_candi) == 1) {
              c2p_prob[1] <- 1 / 3 * swap_two_prob
            }
            p2c_prob[2] <- 1 / 10 # Q(R_c|T_c)
            c2p_prob[2] <- prop_trans_func_prob # Q(R_p|T_p)
            p2c_prob[3] <- 1 # Q(\theta_c|T_c, R_c) no root node is introduced.
            c2p_prob[3] <- 1 # Q(\theta_p|T_p, R_p)

            nume <- sum((prop_sample_prob)) + sum(log(p2c_prob))
            deno <- sum((curr_sample_prob)) + sum(log(c2p_prob))
            acce_prob <- exp(nume - deno)
            ratio <- runif(1)
            if (ratio <= acce_prob) {
              n <- n + 1
            }
          }
          #############################################################
          #  additional move 1: reorder the input and output variables  for XOR relation
          if (uu >= 4 / 7 && uu < 5 / 7 && current_trans_func_matrix[update_order[k], parent_of_update[1]] == 9) {
            prop_incid_matrix <- current_incid_matrix
            sample.child <- parent_of_update[sample.int(2, 1)]
            for (i in 1:2) {
              prop_incid_matrix[update_order[k], parent_of_update[i]] <- 0
            }
            parent_parent <- numeric()
            mm <- 0
            for (i in 1:num.node) {
              if (current_incid_matrix[sample.child, i] > 0) {
                mm <- mm + 1
                parent_parent[mm] <- i
              }
            }
            if (length(parent_parent) == 0) # selected child should have no parents before this operation
              {
                prop_incid_matrix[sample.child, update_order[k]] <- 1
                prop_incid_matrix[sample.child, setdiff(parent_of_update, sample.child)] <- 1
                prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
                if (check_ances_matrix(prop_ances_matrix) == 0) # make sure no directed cycles
                  {
                    prop_trans_func_matrix <- current_trans_func_matrix
                    func_order <- current_trans_func_matrix[update_order[k], parent_of_update[1]]
                    for (i in 1:2) {
                      prop_trans_func_matrix[update_order[k], parent_of_update[i]] <- 0
                    }
                    prop_trans_func_matrix[sample.child, update_order[k]] <- func_order
                    prop_trans_func_matrix[sample.child, setdiff(parent_of_update, sample.child)] <- func_order
                    xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix)
                    prop_sample <- xxx[[2]]
                    prop_post <- xxx[[1]][length(xxx[[1]])]

                    nume <- sum((prop_post))
                    deno <- sum((current_post))
                    acce_prob <- exp(nume - deno)
                    ratio <- runif(1)
                    if (ratio <= acce_prob) {
                      n <- n + 1
                    }
                  }
              }
          }
          ###################################################
          # additional move 2: reverse one arc among gene triplet
          if (uu >= 5 / 7 && uu < 6 / 7) {
            reverse_pare <- parent_of_update[sample.int(2, 1)]
            remain_pare <- setdiff(parent_of_update, reverse_pare)
            ###################
            reverse_node_pare <- numeric()
            j <- 1 # find the already existing parent for reverse_pare
            for (i in 1:ncol(current_incid_matrix)) {
              if (current_incid_matrix[reverse_pare, i] != 0) {
                reverse_node_pare[j] <- i
                j <- j + 1
              }
            }
            if (length(reverse_node_pare) == 0) # currently no parent
              {
                prop_incid_matrix <- current_incid_matrix
                prop_incid_matrix[update_order[k], reverse_pare] <- 0
                prop_incid_matrix[reverse_pare, update_order[k]] <- 1
                prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
                if (check_ances_matrix(prop_ances_matrix) == 0) # make sure no loops
                  {
                    prop_trans_func_matrix <- current_trans_func_matrix
                    prop_trans_func_matrix[update_order[k], reverse_pare] <- 0
                    ###############################
                    prop_trans_func_matrix[update_order[k], remain_pare] <- 10 + sample.int(2, 1) # try no use prior
                    prop_trans_func_matrix[reverse_pare, update_order[k]] <- 10 + sample.int(2, 1)
                    xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix)
                    prop_sample <- xxx[[2]]
                    prop_post <- xxx[[1]][length(xxx[[1]])]

                    nume <- sum((prop_post))
                    deno <- sum((current_post))
                    acce_prob <- exp(nume - deno)
                    ratio <- runif(1)
                    if (ratio <= acce_prob) {
                      n <- n + 1
                    }
                  }
              }
            if (length(reverse_node_pare) == 1) # currently one parent
              {
                prop_incid_matrix <- current_incid_matrix
                prop_incid_matrix[update_order[k], reverse_pare] <- 0
                prop_incid_matrix[reverse_pare, update_order[k]] <- 1
                prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
                if (check_ances_matrix(prop_ances_matrix) == 0) # make sure no loops
                  {
                    prop_trans_func_matrix <- current_trans_func_matrix
                    prop_trans_func_matrix[update_order[k], reverse_pare] <- 0
                    ###############################
                    prop_trans_func_matrix[update_order[k], remain_pare] <- 10 + sample.int(2, 1) # no use prior
                    func_order <- sample.int(10, 1)
                    prop_trans_func_matrix[reverse_pare, update_order[k]] <- func_order
                    prop_trans_func_matrix[reverse_pare, reverse_node_pare] <- func_order
                    xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix)
                    prop_sample <- xxx[[2]]
                    prop_post <- xxx[[1]][length(xxx[[1]])]

                    nume <- sum((prop_post))
                    deno <- sum((current_post))
                    acce_prob <- exp(nume - deno)
                    ratio <- runif(1)
                    if (ratio <= acce_prob) {
                      n <- n + 1
                    }
                  }
              }
          }
          ########################################
          # additional move 3: only change boolean functions among gene triplets  note: the network may remain the same
          if (uu >= 6 / 7) {
            prop.prob <- runif(1)
            aa <- 0
            if (prop.prob >= prop.ratio) {
              if (nrow(triplet.prior.set) >= 1 && ncol(triplet.prior.set) > 1) {
                candidate.pare.set <- list()
                for (jj in 1:nrow(triplet.prior.set)) {
                  if (length(intersect(parent_of_update, triplet.prior.pare[jj, ])) == 2) {
                    aa <- aa + 1
                    candidate.pare.set[[aa]] <- triplet.prior.set[jj, ]
                  }
                }
                if (length(candidate.pare.set) > 0) {
                  score <- numeric()
                  for (jj in 1:length(candidate.pare.set)) {
                    score[jj] <- candidate.pare.set[[jj]][5]
                  }
                  bbb <- matrix(nrow = length(candidate.pare.set), ncol = 6)
                  for (jjj in 1:length(candidate.pare.set)) {
                    bbb[jjj, ] <- candidate.pare.set[[jjj]]
                  }

                  max.score.candidate <- bbb[bbb[, 5] == max(score), ]
                  if (is.numeric(max.score.candidate) == T) {
                    add_parent <- max.score.candidate
                  }
                  if (is.matrix(max.score.candidate) == T) {
                    if (nrow(max.score.candidate) > 1) {
                      add_parent <- max.score.candidate[sample.int(nrow(max.score.candidate), 1), ]
                    }
                  }
                }
              }
            }
            prop_incid_matrix <- current_incid_matrix
            prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
            prop_trans_func_matrix <- current_trans_func_matrix
            if (aa > 0) {
              for (i in 1:2) {
                prop_trans_func_matrix[update_order[k], parent_of_update[i]] <- add_parent[4]
              }
            }
            if (prop.prob < prop.ratio || aa == 0) {
              func_order <- sample.int(10, 1)
              for (i in 1:2) {
                prop_trans_func_matrix[update_order[k], parent_of_update[i]] <- func_order
              }
            }
            xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix)
            prop_sample <- xxx[[2]]
            prop_post <- xxx[[1]][length(xxx[[1]])]
            nume <- sum(prop_post)
            deno <- sum(current_post)
            acce_prob <- exp(nume - deno)
            ratio <- runif(1)
            if (acce_prob != 1 && ratio <= acce_prob) {
              n <- n + 1
            }
          }
        }
        #######################################
        #######################################
        if (old == n) # no update
          {
            iter <- iter + 1
            xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = current_trans_func_matrix)
            current_sample <- xxx[[2]]
            Sample_Matrix[iter, ] <- current_sample
            all_logpost[iter] <- xxx[[1]][length(xxx[[1]])]
            Incidence_Matrix[[iter]] <- current_incid_matrix
            Ancestor_Matrix[[iter]] <- current_ances_matrix
            Trans_Func_Matrix[[iter]] <- current_trans_func_matrix
          }
        if (old != n) # when either T or F changes,  recount the number
          {
            iter <- iter + 1
            jump_point[n] <- iter

            logpost[n] <- sum(prop_post)
            Sample_Matrix[iter, ] <- prop_sample
            all_logpost[iter] <- sum(prop_post)
            Incidence_Matrix[[iter]] <- prop_incid_matrix
            Ancestor_Matrix[[iter]] <- prop_ances_matrix
            Trans_Func_Matrix[[iter]] <- prop_trans_func_matrix
          }
      } # end of updating
    } # end of num_update
    ##########################
    # num<-numeric()  # summarize the result
    # for (i in 2:length(jump_point))
    #     num[i-1]<-jump_point[i]-jump_point[i-1]
    # num[length(jump_point)]<-iter+1-jump_point[length(jump_point)]
    # for (i in 1:length(num))
    #  if (num[i]==max(num))
    #     max_num<-i
    ######################################
    # true_root<-numeric(); j<-0
    # for ( i in 1:nrow(true_network))
    #   if (sum(true_network[i,])==0)
    #   {
    #     j<-j+1; true_root[j]<-i
    #   }
    # common_links<-list(); common_root<-list(); model_root<-list()
    # non_inferred_links<-numeric(); distance_matrix1<-list(); false_links<-numeric()
    # correct_links_ratio1<-numeric(); correct_links_ratio2<-numeric()
    # for (i in 1:length(Incidence_Matrix))
    #   {
    #   k<-0 ; root_node<-numeric()
    #   for (j in 1:nrow(Incidence_Matrix[[i]]))
    #     if (sum(Incidence_Matrix[[i]][j,])==0)
    #       {
    #         k<-k+1; root_node[k]<-j
    #       }
    #   model_root[[i]]<-root_node
    #   common_root[[i]]<-intersect(true_root, root_node)
    #   distance_matrix1[[i]]<-true_incid_matrix-Incidence_Matrix[[i]]
    #   bb<-distance_matrix1[[i]]
    #   cc<-bb[bb>0]; dd<-bb[bb<0]
    #   non_inferred_links[i]<-sum(cc)
    #   false_links[i]<-abs(sum(dd))
    #   correct_links_ratio1[i]<-(sum(Incidence_Matrix[[i]])-false_links[i]+length(common_root[[i]]))/(sum(Incidence_Matrix[[i]])+length(root_node))  # correct_links_ratio=correctly inferred links/inferred total links
    #   correct_links_ratio2[i]<-(sum(true_incid_matrix)-non_inferred_links[i]+length(common_root[[i]]))/(sum(true_incid_matrix)+length(true_root))
    #   }
    # ###########################
    # total_func<-0        # report inferred structure and functions
    # for (i in 1:nrow(true_incid_matrix))
    #   if (sum(true_incid_matrix[i,])>0)
    #     total_func<-total_func+1
    # non_inferred_SF<-numeric(); SF_ratio<-numeric();  distance_matrix2<-list()
    # for (i in 1:length(Trans_Func_Matrix))
    #   {
    #   distance_matrix2[[i]]<-true_network-Trans_Func_Matrix[[i]]
    #   tt<-0
    #   for (j in 1:nrow(distance_matrix2[[i]]))
    #     {
    #       bb<-distance_matrix2[[i]][j,]
    #       if(max(bb)>0)
    #         tt<-tt+1
    #     }
    #   non_inferred_SF[i]<-tt
    #   SF_ratio[i]<-(total_func-tt+length(common_root[[i]]))/(total_func+length(true_root))
    #   }
    #  All_Trans_Func_Matrix[[iii]]=Trans_Func_Matrix
    #  All_Logpost[[iii]]=all_logpost
    #  All_Correct_Rate[[iii]]=SF_ratio
  } # end of iii
  return(list(
    networks = Trans_Func_Matrix,
    log_posterior = all_logpost
    # correct_rate = correct_links_ratio1
  ))
}
