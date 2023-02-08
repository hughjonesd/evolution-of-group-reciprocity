

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
library(purrr)

#' Run a single simulation with specific parameters
#'
#' @param b benefit of being cooperated with
#' @param c cost of cooperation
#' @param G Group size
#' @param n_groups number of groups
#' @param T_rounds number of rounds per generation
#' @param initial_thresholds initial vector of thresholds k, size G*n_groups
#'   A `k >= 1` is a selfish type
#' @param k initial threshold of GR types. This overrides `initial_thresholds` 
#'   and gives a vector of `prop_gr` `k`, `1-prop_gr` `1` (i.e. selfish types),
#'   randomized across groups
#' @param prop_gr initial proportion of GR types.
#' @param q what proportion of the new generation reproduces by "local" fitness
#'   i.e. from the local group proportionally to fitness. Note that groups are
#'   always randomized between generations.
#' @param generations how many generations to run for (max; if we hit fixation
#'   we stop early with a warning).
#'
#' @return Vector of thresholds k after all generations, or at fixation
run_simulation <- function (
                            b, 
                            c, 
                            G, 
                            n_groups, 
                            T_rounds, 
                            initial_thresholds = NULL, 
                            k = NULL, 
                            prop_gr = 0.5,
                            q = 0, 
                            generations = 100
                            ) {
  stopifnot(G > 1, n_groups > 1, T_rounds > 1, c >=0, b >= c, q >= 0, q <= 1,
            generations >= 1)
  if (! is.null(initial_thresholds)) {
    stopifnot(length(initial_thresholds) == G * n_groups,
            all(initial_thresholds >= 0))
  } else {
    stopifnot(k >= 0, k <= 1, prop_gr >= 0, prop_gr <= 1)
  }
  
  params <- list(b = b, c = c, G = G, n_groups = n_groups, T_rounds = T_rounds,
                 q = q)
  thresholds <- if (is.null(k)) {
                  initial_thresholds
                } else {
                  thr <- rep(1, G*n_groups)
                  thr[1:ceiling(prop_gr * G * n_groups)] <- k
                  sample(thr)
                }
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
  unchanged <- FALSE
  for (t in 1:T_rounds) {
    # if cooperation is unchanged from last time, result will be unchanged too:
    if (! unchanged) {
      result <- run_one_round(coop_last_round = coop, thresholds = thresholds, 
                            params = params)
    }
    payoffs <- payoffs + result$payoffs
    unchanged <- all(result$coop == coop)
    coop <- result$coop
  }
  
  return(payoffs)
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

  # calculate past average cooperation of groups against every indiv
  group_list <- lapply(seq(1, by = G, length = n_groups),
                         function (ix) coop_last_round[ix:(ix+G-1), ]
                       )
  mean_group_coop <- vapply(group_list, colMeans, FUN.VALUE = numeric(N))
  mean_group_coop <- t(mean_group_coop)
  
  stopifnot( dim(mean_group_coop) == c(n_groups, N) )

  # note that the strict inequality means a threshold of 1 is a selfish type
  # you can use a threshold of e.g. 0.99 for a true 1-threshold GR
  coop_with_group <- apply(mean_group_coop, 1, function (mgc) mgc > thresholds)
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


#' Run the selection process after a generation
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
  if (q > 0) { # speed optimization
    dfr <- data.frame(groups = groups, thresholds = thresholds, payoffs = payoffs)
    dfr <- dfr |> 
           group_by(groups) |> 
           mutate(
             wg_thresholds = sample(thresholds, replace = TRUE, prob = payoffs)
           )
    # ... but only apply it for the within-group people
    new_thresholds[within_group] <- dfr$wg_thresholds[within_group]
  }
  
  new_thresholds[!within_group] <- sample(thresholds, 
                                          size = sum(!within_group), 
                                          replace = TRUE, 
                                          prob = payoffs)
  
  # always randomize groups
  new_thresholds <- sample(new_thresholds)
  
  return(new_thresholds)
}


stop("Tests below")


# experiments to run:
# * vary n_groups
# * vary T_rounds
# * vary q
# 

params <- data.frame(b = 5, c = 1, G = 8, n_groups = 50, k = 0.5, T_rounds = 20,
                     q = 0, generations = 100)

params <- bind_rows(
  Basic    = params,
  `High k = 0.7`    = params |> mutate(k = 0.7),
  `High G = 16`     = params |> mutate(G = 16),
  `Few groups = 20` = params |> mutate(n_groups = 20),
  `Short T = 10`     = params |> mutate(T_rounds = 10),
  `Local selection q = 0.25` = params |> mutate(q = 0.25),
  .id = "variant"
)

n_reps <- 30
params <- params[rep(1:nrow(params), n_reps), ]

prop_gr <- params |> 
           select(-variant) |> 
           purrr::pmap(run_simulation) |> 
           purrr::map_dbl(\(x) mean(x < 1))

params$prop_gr <- prop_gr
params$rep <- 1:n_reps

library(ggplot2)

ggplot(params, aes(variant, prop_gr, color = variant)) + 
  geom_point(alpha = 0.5) + 
  stat_summary(fun = mean, shape = "cross", size = 1.3) + 
  labs(y = "Prop. GR types") +
  theme_minimal() +
  theme(legend.position = "none", text = element_text(size = 12))

# a simulation with large groups
# starting with high k and many initial GR types helps, since then
# there are likely to be some supraliminal groups, but not many
# free riders in them
big <- run_simulation(b = 5, c = 1, 
                      k = 0.9, prop_gr = 0.9,
                      G = 40, n_groups = 20, 
                      T_rounds = 20,  q = 0, generations = 1000)


# a simulation with many starting thresholds and 25% selfish
# this typically reduces to one or two types, which suggests there's a coordination
# advantage to being the same type as others
poly <- run_simulation(b = 5, c = 1, 
                      initial_thresholds = sample(c(runif(300), rep(1, 100))),
                      G = 8, n_groups = 50, 
                      T_rounds = 20,  q = 0, generations = 1000)

polys <- replicate(20, {
          run_simulation(b = 5, c = 1, 
                        initial_thresholds = sample(c(runif(300), rep(1, 100))),
                        G = 8, n_groups = 50, 
                        T_rounds = 20,  q = 0, generations = 500)
        })
apply(polys, 2, \(x) length(unique(x)))
min_and_max <- apply(polys, 2, range)
data.frame(min = min_and_max[1,], max = min_and_max[2,]) |> 
  arrange(max, min) |> 
  mutate(rep = 1:20) |> 
  ggplot() + 
    geom_segment(aes(x = rep, xend = rep, y = min, yend = max)) +
    geom_point(aes(x = rep, y = min)) + 
    geom_point(aes(x = rep, y = max)) + 
    labs(y = "Thresholds", x = "Simulation")



# approximating the model: many groups, small groups (for speed), largish T,
# vary k around c/b, many generations

model <- data.frame(
  k = rep(c(1/4, 1/3, 1/2) + 0.01, each = 3),
  b = rep(2:4, 3),
  c = 1,
  T_rounds = 50,
  generations = 1000,
  n_groups = 150,
  G = 12
)
model$prop_gr <- model$k

prop_gr <- model |> 
           purrr::pmap(run_simulation) |> 
           purrr::map_dbl(\(x) mean(x < 1))
