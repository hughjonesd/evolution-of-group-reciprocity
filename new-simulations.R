

# simulations for evolution of group reciprocity
# each generation
# Groups of size G, each total N = nG
# T_rounds rounds, everyone plays a PD against everyone
# cost of cooperation -c, benefit of being cooperated with b
#
# Selfish types always defect
# GR types cooperate in round 1, then cooperate with group g
# if the proportion of people who cooperated with you from group g
# in round t - 1 is above a threshold k \in (0,1).
# 
# Between generations, reproduction is according to fitness
# Assumption 1: varying within & between-group fitness
# proportion q are chosen from within-group types (proportionally to fitness)
# proportion 1-q from the whole population (proportionally to fitness)
# q = 0 is like the model; with q = 1, selfish types always win
# 
# 
# Other things to explore: 
#  - 3 thresholds, including k > 1 for a selfish type?
#  
#  

library(dplyr)

#' Title
#'
#' @param G Group size
#' @param b benefit of being cooperated with
#' @param c cost of cooperation
#' @param n_groups number of groups
#' @param T_rounds number of rounds per generation
#' @param initial_thresholds initial vector of thresholds k, size G*n_groups
#'
#' @return
#' @export
#'
#' @examples
run_simulation <- function (G, b, c, n_groups, T_rounds, initial_thresholds, q,
                            generations = 1000L) {
  stopifnot(G > 1, n_groups > 1, T_rounds > 1, c >=0, b >= c, q >= 0, q <= 1,
            length(initial_thresholds) == G * n_groups,
            all(initial_thresholds >= 0), 
            generations >= 1)
  
  params <- list(b = b, c = c, G = G, n_groups = n_groups, T_rounds = T_rounds,
                 q = q)
  
  thresholds <- initial_thresholds
  for (gen in 1:generations) {
    payoffs <- run_one_generation(thresholds = thresholds, params = params)
    thresholds <- run_selection_process(thresholds = thresholds, 
                                        payoffs = payoffs, params = params)
    
    if (length(unique(thresholds)) == 1L) {
      warning("Reached fixation at generation ", gen)
      return(thresholds)
    }
  }
  
  return(thresholds)
}

#' Run a single round of cooperation
#'
#' @param coop_last_round N*N matrix of cooperation (T) or 
#'   defection (F) in previous round. Row i, j represents whether i cooperated
#'   with j.
#' @param thresholds G*n_groups vector of thresholds
#' @param params List of simulation global parameters: b, c, n_groups, G, T_rounds
#'
#' @return list containing `coop` N*N matrix of cooperation or defection this round;
#'   and `payoffs` N-vector of this-round payoffs
run_one_round <- function(coop_last_round, thresholds, params) {
  # Rows 1-G are the behaviour of the first group
  # Rows G+1 - 2G are the behaviour of the second group etc.
  G <- params$G
  n_groups <- params$n_groups
  N <- G * n_groups
  group_index <- rep(1:n_groups, each = G)

  # returns a n_groups length vector of the proportion who cooperated
  # in each group, given a N-vector of individual cooperation
  calculate_mean_group_coop <- function (coop_vector) {
    tapply(coop_vector, group_index, mean)
  }
  
  # calculate cooperation against every individual
  mean_group_coop <- apply(coop_last_round, 2, calculate_mean_group_coop)
  stopifnot( dim(mean_group_coop) == c(n_groups, N) )

  coop_with_group <- apply(mean_group_coop, 1, function (mgc) mgc >= thresholds)
  # row i, column j says whether individual i should cooperate with
  # everyone in group j

  # to create the full NxN matrix, we just repeat the group columns G times each:
  coop <- coop_with_group[, group_index]
  
  # each individual loses c from cooperating and gains b from being coop-ed with
  # we normalize so that the sucker's payoff is zero
  payoffs <- N * params$c - 
             rowSums(coop) * params$c +
             colSums(coop) * params$b 
  
  list(payoffs = payoffs, coop = coop)
}

#' Run a generation of T_rounds rounds
#'
#' @param thresholds,params As in `run_one_round`
#'
#' @return A N-vector of total payoffs from the rounds
run_one_generation <- function (thresholds, params) {
  T_rounds <- params$T_rounds
  N <- params$n_groups * params$G
  
  payoffs <- rep(0, N)
  # we start the process by creating a fake round zero
  # where everybody cooperated. Then only "selfish" types with
  # a threshold above 1 will not cooperate in round 1
  coop <- matrix(TRUE, N, N)
  for (t in 1:T_rounds) {
    result <- run_one_round(coop_last_round = coop, thresholds = thresholds, 
                            params = params)
    payoffs <- payoffs + result$payoffs
    coop <- result$coop
  }
  
  return(payoffs)
}


#' Run the selection process
#'
#' @param thresholds Existing thresholds N-vector
#' @param payoffs Payoffs from the T rounds
#' @param params As above. In particular, `q` is important (prob of within-group
#'   selection)
#'
#' @return New thresholds N-vector
run_selection_process <- function (thresholds, payoffs, params) {
  q <- params$q
  # with probability q, you sample within the group, weighted by payoffs
  # otherwise you sample within the whole population, weighted by payoffs
  
  n_groups <- params$n_groups
  G <- params$G
  N <- n_groups * G
  groups <- rep(1:n_groups, each = G)
  
  within_group <- runif(N) <= q
  new_thresholds <- numeric(N)
  
  # calculate within-group selection for everyone...
  dfr <- data.frame(groups = groups, thresholds = thresholds, payoffs = payoffs)
  dfr <- dfr |> 
         group_by(groups) |> 
         mutate(
           wg_thresholds = sample(thresholds, replace = TRUE, prob = payoffs)
         )
  # ... but only apply it for the within-group people
  new_thresholds[within_group] <- dfr$wg_thresholds[within_group]
    
  new_thresholds[!within_group] <- sample(thresholds, 
                                          size = sum(!within_group), 
                                          replace = TRUE, 
                                          prob = payoffs)
  
  return(new_thresholds)
}


params <- list(n_groups = 10, G = 20, b = 2, c = 1, T_rounds = 20, q = 0, 
               generations = 1000)
thresholds <- rep(c(0.8,  1.01), each = params$G * params$n_groups / 2)
thresholds <- sample(thresholds)

# run_one_round(coop_last_round = matrix(TRUE, 9, 9), thresholds = thresholds, params = params)
# payoffs <- run_one_generation(thresholds = thresholds, params = params)
# new_thresholds <- run_selection_process(thresholds, payoffs, params)

result <- run_simulation(G = 20, b = 3, c = 1, n_groups = 10, 
               initial_thresholds = thresholds, 
               T_rounds = 20, 
               q = 0, 
               generations = 50)
