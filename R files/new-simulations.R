

# simulations for evolution of group reciprocity
# each generation
# n_groups groups of size G, total N = n_groups * G
# T_rounds rounds, everyone plays a PD against everyone
# cost of cooperation c, benefit of being cooperated with b
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

## setup ====
library(dplyr)
library(purrr)
library(furrr)
library(ggplot2)
library(glue)

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
  changed <- TRUE
  for (t in 1:T_rounds) {
    # if cooperation is unchanged from last time, result will be unchanged too:
    if (changed) {
      result <- run_one_round(coop_last_round = coop, thresholds = thresholds, 
                            params = params)
    } else {
      # we can just add the payoffs times number of remaining rounds, then break
      payoffs <- payoffs + result$payoffs*(T_rounds + 1 - t)
      break
    }
    payoffs <- payoffs + result$payoffs
    # max() appears faster than any()
    changed <- max(result$coop != coop) 
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

# The idea: rather than running multiple generations,
# we just start with a given prop_gr and run 1 generation
# many times. This gives us the state space. In particular if
# we take the mean prop GR after 1 generation then we could
# think of this as a large metapopulation of many populations
# (who would all reform after 1 generation); and draw the phase
# diagram so that prop GR is increasing when the mean final prop is
# above the init prop
# 

#' Map the state space for given parameters over 1 generation
#'
#' @param ... Parameters passed to [run_simulation()]. You can include
#'   multiple values of parameters, they'll be crossed using [expand.grid()].
#' @param nreps How many repetitions of each parameter combination to run?
#' @param prop_gr Initial proportions of group reciprocators. Default
#'   is 0.05, 0.1, 0.15, ... 0.95
#'
#' @return A data frame giving mean and s.e. of the proportion of group
#'   reciprocators after 1 generation, starting from different initial
#'   proportions,
#'   for each parameter combination. It also gives the mean increment, i.e.
#'   the final proportion of group reciprocators minus the initial proportion.
map_state_space <- function (..., nreps = 20, prop_gr = seq(0.05, 0.95, 0.05)) {
  params <- expand.grid(..., 
                        generations = 1,
                        prop_gr = prop_gr, 
                        experiment = seq_len(nreps))
  params$prop_gr_final <- params |> 
                          select(-experiment) |> 
                          purrr::pmap(run_simulation, .progress = TRUE) |> 
                          purrr::map_dbl(\(x) mean(x < 1))
  
  params |> 
    group_by(across(c(-experiment, -prop_gr_final, -generations))) |> 
    summarize(
      mean_gr_final = mean(prop_gr_final),
      incr_gr_final = mean_gr_final - prop_gr,
      se_gr_final = sd(prop_gr_final)/sqrt(nreps)
    )
}


stop("Tests below")

## Basic tests ====

# Idea: check that having finite number of groups, and finite T, still allows
# GR to evolve; and that the basic k > c/b "comparative statics" still hold.
# (I.e. we'd expect higher k, and lower c/b ratio, to make evolution of GR
# more likely.)
# 
# Plan: run X experiments; in each run 100 generations
#  take the proportion of Gr at the end; 

set.seed(27101975)

n_simulations <- 100
group_size_G <- 8
n_generations <- 500

basic_params <- expand.grid(
                  b = c(2, 5), 
                  c = 1, 
                  G = group_size_G, 
                  n_groups = c(10, 50, 200),
                  T_rounds = c(10, 20, 50), 
                  k = c(0.4, 0.8), 
                  generations = n_generations,
                  experiment = 1:n_simulations 
                ) |> 
  filter(
    (n_groups < 200 & T_rounds < 50) | (n_groups == 200 & T_rounds == 50) 
  )


basic_params$prop_gr_final <- basic_params |> 
                        select(-experiment) |> 
                        purrr::pmap(run_simulation, .progress = TRUE) |> 
                        purrr::map_dbl(\(x) mean(x < 1))

today <- Sys.Date()
saveRDS(basic_params, file = glue::glue("R files/basic-params-{today}.Rds"))

basic_results <- basic_params |> 
                 group_by(b, n_groups, T_rounds, k) |>
                 summarize(
                   mean_gr = mean(prop_gr_final),
                   prop_fix_0 = mean(prop_gr_final == 0),
                   prop_fix_1 = mean(prop_gr_final >= 1),
                   prop_fix = prop_fix_0 + prop_fix_1
                 )

basic_results |> 
  mutate(
    `N groups` = factor(n_groups),
    `T` = factor(T_rounds),
    Condition = glue::glue("{T_rounds} rounds, {n_groups} groups"),
    `c/b` = factor(1/b, levels = c(0.2, 0.5)),
    k = factor(k)
  ) |> 
  ggplot(aes(x = `c/b`, y = mean_gr, shape = k, color = k)) + 
    facet_grid(vars(`N groups`), vars(`T`), labeller = label_both) + 
    geom_point(size = 3, position = position_dodge(width = 0.1)) + 
    theme_linedraw() + 
    labs(
      y = "Prop. GR types", 
      subtitle = glue::glue("{n_simulations} simulations, G = {group_size_G}, max {n_generations} generations")
    )
    

ggsave("R files/basic-experiment.jpeg", width = 5, height = 4)

# TODO:
# - regenerate the above figure with more generations & replications;
# - add a row with n_groups = 200, a column with T = 50
# - check we approach the theoretical result?
# - simulations with a power function proportional response:
#   you cooperate with probability x^alpha where x is the prop in the group
#   that cooperated with you; alpha > 0.
# - theory of power functions. Does it collapse to threshold k = 1 when
#   alpha >= 1? more generally, for any line below the diagonal?
# - simulate evolution of thresholds?
#   - either mutation during the sim, or simulate invasion by a mutant
#     starting from an eqm distribution of types.


## Evolution from many thresholds ====
 
set.seed(27101975)
plan(multisession, workers = 6)

n_simulations <- 100
group_size_G <- 50
n_generations <- 2000
n_groups <- 20
T_rounds <- 20

fixation_params <- expand.grid(
  b = c(2, 5),
  c = 1,
  G = group_size_G,
  n_groups = n_groups,
  T_rounds = T_rounds,
  generations = n_generations,
  experiment = 1:n_simulations
) |> 
  rowwise() |> 
  mutate(
    initial_thresholds = list(
      c(runif(n_groups * G * 0.9), 
        rep(1, n_groups * G * 0.1))
    )
  ) |> 
  ungroup()

fixation_params$result <- fixation_params |> 
                        select(-experiment) |> 
                        furrr::future_pmap(run_simulation, .progress = TRUE,
                                           .options = furrr::furrr_options(seed = TRUE))
                        

fixation_results <- fixation_params |> 
  rowwise() |> 
  mutate( 
    result_n = n_distinct(result),
    result_min = min(result),
    result_max = max(result),
    result_mean = mean(result)
  ) |> 
  ungroup()

fixation_results |> 
  mutate(
    `Distinct thresholds` = factor(result_n),
    `Mean threshold` = result_mean,
    `b/c` = factor(b)
  ) |> 
  ggplot(aes(x = `b/c`, y = `Mean threshold`)) + 
    geom_violin(draw_quantiles = 1:3/4) +
    geom_point(aes(color = `Distinct thresholds`), alpha = 0.5) +
    theme_linedraw() +
    labs(
      subtitle = glue::glue("{n_simulations} simulations, max {n_generations} generations, 
                            {n_groups} groups, {T_rounds} rounds")
    )


## Checking for new entrants, starting with selfish ====
## 

set.seed(27101975)
plan(multisession, workers = 6)

n_simulations <- 100
n_generations <- 1000
T_rounds <- 50

entrant_params <- expand.grid(
  b = 5,
  c = 1,
  G = c(10, 20),
  n_groups = c(20, 40),
  T_rounds = T_rounds,
  generations = n_generations,
  experiment = 1:n_simulations,
  prop_gr = 0.2,
  k = seq(0.1, 0.9, 0.1)
)  

entrant_params$result <- entrant_params |> 
                         select(-experiment) |> 
                         furrr::future_pmap(run_simulation, .progress = TRUE,
                                           .options = furrr::furrr_options(seed = TRUE))
 
entrant_results <- entrant_params |> 
  rowwise() |> 
  mutate( 
    result_n = n_distinct(result),
    result_min = min(result),
    result_max = max(result),
    result_mean = mean(result)
  ) |> 
  ungroup()

entrant_results |> 
  summarize(.by = c(G, n_groups, k),
    prop_survived = mean(result_mean < 1),
    prop_fixated = mean(result_n == 1 & result_mean < 1)
  ) |> 
  ggplot(aes(k, prop_survived)) +
    geom_col() +
    geom_col(aes(y = prop_fixated), fill = alpha("orange", 0.5)) +
    facet_grid(vars(G), vars(n_groups), labeller = label_both)
  

## Finding the state space ==== 

# This is old, not sure if relevant


ss <- map_state_space(b = c(2, 5, 8), c = 1, k = c(0.4, 0.6, 0.8), G = 8, 
                      n_groups = 20, T_rounds = c(10, 20), 
                      nreps = 100, )


ss |> 
  mutate(
    inc_gr = mean_gr_final - prop_gr
  ) |> 
  ggplot(aes(x = prop_gr, y = inc_gr, color = factor(k), group =k)) + 
    geom_abline(intercept = 0, slope = 0) +
    geom_line() +
    geom_pointrange(aes(ymax = inc_gr + 1.96*se_gr_final, 
                        ymin = inc_gr - 1.96*se_gr_final), size = 0.3) +
    facet_grid(vars(T_rounds), vars(b), labeller = label_both) +
    theme_linedraw()
  

ss2 <- map_state_space(b = 8, c = 1, k = c(0.4, 0.8), G = c(8, 16, 24), 
                      n_groups = c(10, 20), T_rounds = 10, nreps = 100)

ss2 |> 
  mutate(
    inc_gr = mean_gr_final - prop_gr
  ) |> 
  ggplot(aes(x = prop_gr, y = inc_gr, color = factor(k), group = k)) + 
    geom_abline(intercept = 0, slope = 0) +
    geom_line() +
    geom_pointrange(aes(ymax = inc_gr + 1.96*se_gr_final, 
                        ymin = inc_gr - 1.96*se_gr_final), size = 0.3) +
    facet_grid(vars(n_groups), vars(G), labeller = label_both) +
    theme_linedraw()
