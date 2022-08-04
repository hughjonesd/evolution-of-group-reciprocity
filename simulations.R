
# each generation has a strategy and a group ID (G groups of N/G indivs each, say N = 3*4*5 to allow G = 2,3,4,5)
# for T rounds, 
  # indivs are paired (more likely within group), and:
  # indivs play cooperate or defect (payoffs t > 1 > 0 > p) according to their strategies:
  # Selfish: always defects
  # Altruist: always coops
  # GEneralized Reciprocity: plays same as own last "experience" (what someone played against him); starts with coop
  # Group Reciprocity: plays same as own last "experience" from current partner's group; cooperates if no experience
  # Group-Group reciprocity: plays same as own group's last "experience" from current partner's group (randomize if nec)

# sum payoffs
# replicator dynamics within group: with some prob, indivs change strategy and are more likely to pick most successful


# what we think is happening:
# fitness is a function of group reputation. Groups benefit from a good group reputation. Indivs within a group
# benefit from cheating. So, there is a standard Simpson's paradox / George Price thing happening.
# Predictions:
# selfish types are better off within each group; but groups with more GRs are better off.
# GRs do better if there is more type heterogeneity between, rather than within groups
# between group type heterogeneity should also correlate with *fitness* variance being high between groups (at last period)
#   and with reputation variance between groups also
# groups with many GRs are helped more, and have better average "reputation"/standing

# comparative statics:
# - larger groups mean selfish types win (less heterogeneity/granularity and also public good effect (maybe?))
# - more periods better for GRs (they spend less time being naive and grow a better rep)
# - benefit cost ratio should help GRs, I guess
# - initial type distribution, how much between-group variance

# PLAN:
# write tests for strategies doing what you expect and for fitness doing what you expect
# speedup if need be
# run simulations until fixation, many times; varying the parameters group size, b, periods, initial type dist
# record time till fixation

# also: 
# compare group reciprocity with Outgroup Generalized Recip ("GE") 
# - (needs code rewrite! GE just averages across all groups?)
# check what happens with copying/mutations within/between groups (some level of within-ness may help, by
#   increasing the group heterogeneity)


# rough empirics so far:
# - large groups find it much harder to have GR fixation
# - b/c ratio does help a bit, but relationship btw this and gsize isn't linear 
# - you need something like b proportion to gsize^2 ? E.g. gsize 60 b 800 worked
# - with memory of 1, by end of 25 periods, nobody is cooperating; but GRs can still have an advantage sometimes, 
#   presumably in middle periods?
# - memory = 1 harms GRs; also introduces more fluctuations into GR/selfish fitness ratio
# - altruists themselves can survive if b/c is great enough and if there are some GR types
# - and altruists might help GR types survive by creating between group variance
library(tibble)
library(magrittr)
library(testthat)
library(glue)
library(scales)
options(warn = 1)


make_matches <- function (group) {
  ngroups <- length(unique(group))

  # everyone is matched with a random outgroup member.
  pop <- cbind(indiv = 1:length(group), group = group) # matrices are fast
  within_order <- c(replicate(ngroups, sample.int(length(group)/ngroups)))
  pop <- pop[order(within_order), ]
  partner <- pop[1:0 + rep(seq(1, nrow(pop) ,2), each = 2), 1]
  
  partner[order(pop[, 1])]
}



do_report <- function (pop, gen, stats) {
  cat("\n\nGeneration", gen, "\n")
  cat("Share of GRs:", round(mean(pop$strat == 'GR'), 3), "\n")
  group_shares <- tapply(pop$strat == 'GR', pop$group, mean)
  cat("Std dev of groups' shares of GRs: ", round(sd(group_shares), 3), "\n")
  table(pop$strat, Cooperation = pop$coop) %>% prop.table(1) %>% round(2) %>% print
  table(pop$strat, Helped = pop$helped) %>% prop.table(1) %>% round(2) %>% print
  GRSfr <- mean(pop$fitness[pop$strat == 'GR'])/mean(pop$fitness[pop$strat == 'S'])
  cat(glue("GR/selfish fitness ratio: {round(GRSfr, 3)}\n\n"))
  print(stats, digits = 2)
  points(gen, mean(pop$strat == 'S'), col = 'red')
  if (any(pop$strat == 'A')) points(gen, mean(pop$strat == 'A'), col = 'orange', pch = 5)
  # points(gen, mean(pop$strat == 'GE'), col = 'darkgreen', pch = 2)
  points(gen, mean(pop$strat == 'GR'), col = 'blue', pch = 3)
  # points(gen, mean(pop$strat == 'GG'), col = 'navy', pch = 4)
  points(gen, GRSfr/2, col = alpha("grey", 0.3), pch = 19)
}


decide_coop <- function (pop, group_memory, period, memory) {
  coop <- rep(NA, nrow(pop))
  coop[pop$strat == 'S']  <- FALSE
  coop[pop$strat == 'A']  <- TRUE
  if (! any(pop$strat == 'GR')) return(coop)
  mem_len <- min(period, memory)
  n_GR <- sum(pop$strat == 'GR')
  relevant_memories <- matrix(NA, ncol = 3, nrow = n_GR)
  relevant_memories[, 1] <- which(pop$strat == 'GR')
  relevant_memories[, 2] <- pop$partner_group[pop$strat == 'GR']
  relevant_memories <- relevant_memories[rep(seq_len(nrow(relevant_memories)), each = memory), , drop = FALSE]
  relevant_memories[, 3] <- rep(1:memory, n_GR)
  # group_memory has each indiv by row, each partner group by col, 
  # and is last "memory" periods (even when period is less)

  gm <- group_memory[relevant_memories] # always a vector
  gm <- matrix(gm, ncol = memory, byrow = TRUE) 
  # this uses memory limit by N of episodes, rather than simply N last rounds. It is MUCH slower:
  # coop[pop$strat == 'GR'] <- runif(n_GR) < apply(gm, 1, function (x) {
  #   x <- x[!is.na(x)]
  #   x <- x[seq(from = length(x), by = -1, length.out = min(length(x), mem_len))]
  #   mean(x, na.rm = TRUE)
  # })
  # coop[pop$strat == 'GR'] <- runif(n_GR) < rowMeans(gm[, max((period - mem_len), 1):period, drop = FALSE], na.rm = TRUE)
  coop[pop$strat == 'GR'] <- runif(n_GR) < rowMeans(gm, na.rm = TRUE)
  coop[is.na(coop)] <- TRUE
  coop
} 


compute_stats <- function (pop, group_memory, gen) {
  group_memory <- aperm(group_memory, c(1, 3, 2)) # so we can use colMeans
  group_reput <- colMeans(group_memory, dims = 2, na.rm = TRUE)
  is_GR <- pop$strat == "GR"
  prop_GR      <- tapply(is_GR, pop$group, mean)
  group_fitness <- tapply(pop$fitness, interaction(pop$group, pop$strat == "GR"), mean, na.rm = TRUE)
  group_fitness_S <- group_fitness[length(group_fitness)/2 + 1:(length(group_fitness)/2)] 
  group_fitness_GR <- group_fitness[1:(length(group_fitness)/2)] 
  suppressWarnings({
    cor_prop_fit <- cor(prop_GR, group_fitness_GR, use = "complete")
    cor_prop_reput <- cor(prop_GR, group_reput, use = "complete")
    cor_reput_fit <- cor(group_reput, group_fitness_GR, use = "complete")
  })
  between_var <- summary(lm(strat=="GR" ~ factor(group), pop))$r.squared
  
  c(gen, cor_prop_fit, cor_prop_reput, cor_reput_fit, between_var)
}


simulate <- function (
        N                   = 60, 
        G                   = 2, 
        generations         = 100,  # controls the plot size 
        stop_at_generations = FALSE,
        memory              = 5,
        fixation_threshold  = 1,
        periods             = 20, 
        init_probs          = c(S = .5, GR = .5),
        c                   = 1, 
        b                   = 1.5,
        switch_prob         = 0.05,
        report              = 1 # report once per 'report' generations, or FALSE for no reports.
      ) {
  stopifnot(N/G == floor(N/G))
  gsize <- floor(N/G)
  stopifnot(gsize %% 2 == 0) # otherwise we can't match within group
  stopifnot(all(names(init_probs) %in% c('S', 'GE', 'GR', 'GG', 'A')))
  stopifnot(c > 0 && b > c)
  init_strats <- sample(names(init_probs), size = N, replace = TRUE, prob = init_probs)
  
  pop <- tibble(
    group = rep(1:G, each = gsize),
    strat = init_strats,
    memory = rep(NA, N)
  )

  gen <- 1
  # 
  stats <- matrix(NA_real_, nrow = 0, ncol = 5)
  colnames(stats) <- c('gen', 'cor_prop_fit', 'cor_prop_reput', 'cor_reput_fit', 'between_var')
  while (TRUE) {
    if (gen > generations && stop_at_generations) break
    if (gen %% generations == 1) {
      stats <- rbind(stats, matrix(NA_real_, nrow = generations, ncol = ncol(stats)))
      if (is.numeric(report)) {
        plot(1:generations, rep(0, generations), type = 'n', ylim = c(0, 1), ylab = '')
        abline(h = 0.5, lty = 2)
      }
    } 
    group_memory <- array(NA, c(N, G, memory))
    pop$fitness    <- periods * c # so you never get 0
    for (p in 1:periods) {
      pop$matches <- make_matches(pop$group)
      pop$partner_group <- pop$group[pop$matches]
      pop$coop <- decide_coop(pop, group_memory, p, memory)
      pop$helped <- NA
      pop$helped[pop$matches] <- pop$coop
      pop$fitness <- pop$fitness + pop$helped * b - pop$coop * c
      # store experience
      memory_idx <- cbind(pop$matches, pop$group, min(p, memory)) # 3D matrix indexing
      # forget old experience: just rearrange, last row gets overwritten
      if (p > memory && memory > 1) group_memory <- group_memory[, , c(2:memory, 1)] 
      group_memory[ memory_idx ] <- pop$coop
    }
    # record stats
    # generation; reputation of groups (by group_memory); composition of groups; avg fitness of groups*types;
    stats[gen,] <- compute_stats(pop, group_memory, gen)
      
    # and correlate composition/reputation/fitness of groups; type heterogeneity;  
    # at end record type heterogeneity and fitness types
    if (is.numeric(report) && (gen %% generations) %% report == 0) do_report(pop, gen %% generations, stats[gen,])
    
    may_switch <- runif(N) < switch_prob
    pop$strat[may_switch] <- sample(pop$strat, size = sum(may_switch), replace = TRUE, prob = pop$fitness)
    if (max(table(pop$strat))/N >= fixation_threshold) {
      message("Reached fixation threshold")
      tbl <- table(pop$strat)
      avg_stats <- colMeans(stats, na.rm = TRUE)
      return(list(
              winner            = names(tbl)[which.max(tbl)],
              generation        = gen,
              cor_prop_fit      = avg_stats['cor_prop_fit'],
              cor_prop_reput    = avg_stats['cor_prop_reput'],
              cor_reput_fit     = avg_stats['cor_reput_fit'],
              between_var       = avg_stats['between_var']
            ))
    }
    gen <- gen + 1
  }
}



# Simulations -----------------------------------------------------------------------------------------------------

run_sims <- function(
        repetition         = NULL, # for parallel
        nreps              = 10, 
        b                  = c(5, 25), 
        periods            = c(24, 48), 
        memory             = c(12, 24), 
        G                  = c(6, 12), 
        N                  = c(72, 240), 
        fixation_threshold = 0.99
      ) {
  if (is.null(repetition)) repetition = 1:nreps
  params <- tidyr::crossing(
          repetition         = repetition, # vary this first, so we get an early idea of how long things will take
          b                  = b, 
          periods            = periods, 
          memory             = memory, 
          G                  = G, 
          N                  = N, 
          fixation_threshold = fixation_threshold
        )

  params <- params %>% filter(periods >= memory)
  results <- params %>% rowwise %>%  do({
    cat(glue::glue("Running simulation: repetition {.$repetition} of {nreps} 
          b = {.$b}, periods = {.$periods}, memory = {.$periods}, 
          G = {.$G}, N = {.$N}, fixation_threshold = {.$fixation_threshold}
          Starting at {date()}\n\n"))
    as.data.frame(simulate(b = .$b, periods = .$periods, memory = .$memory, G = .$G, N = .$N, 
          generations = 100, report = FALSE, fixation_threshold = .$fixation_threshold))
  })
  bind_cols(params, results)
}

# Tests -----------------------------------------------------------------------------------------------------------


test_decide_coop <- function () {
  test_that("Selfish/altruist doesn't/does cooperate", {
    for (strat in c("S", "A")) {
      pop <- tibble(group = 1, strat = strat, partner_group = 2)
      for (mem in c(NA, 0, 1)) {
        group_memory <- array(mem, c(1, 1, 1))
        expect_identical(decide_coop(pop, group_memory, 1, 5), strat == "A")
      }
    }
  }) 
  test_that("Group reciprocators cooperate based on memory", {
    pop <- tibble(group = 2, strat = "GR", partner_group = 1)
    group_memory <- array(NA, c(1, 1, 5))
    expect_true(decide_coop(pop, group_memory, 1, 5))
    for (length in 1:5) {
      group_memory[1, 1, ] <- c(rep(0, length), rep(NA, 5 - length))
      expect_false(decide_coop(pop, group_memory, length, 5)) 
      group_memory[1, 1, ] <- c(rep(1, length), rep(NA, 5 - length))
      expect_true(decide_coop(pop, group_memory, length, 5))
    }
  })
  for (total in c(4, 6, 8, 25)) for (trues in 1:3) {
    test_that(glue("Group reciprocators cooperate right amount for {trues} in {total} coopn history"), {
      pop <- tibble(group = 2, strat = "GR", partner_group = 1)
      group_memory <- array(NA, c(1, 1, total))
      for (resamples in 1:3) {
        history <- sample(c(rep(1, trues), rep(0, total - trues))) # No point using NAs here
        group_memory[1, 1, ] <- history
        prob_coop <- trues/total
        coop <- replicate(1000, decide_coop(pop, group_memory, total, total))
        expect_lt(mean(coop), qbinom(0.9995, 1000, prob_coop) / 1000)
      }
    })
  }
}


test_make_matches <- function () {
  for (gsize in (1:5)*4) for (ngroups in 2:5) for (tries in 1:10){
    test_that(glue("make_matches works for groupsize {gsize} and ngroups {ngroups}"), {
      groups <- rep(1:ngroups, each = gsize)
      partner <- make_matches(groups)
      expect_true(all(1:length(groups) != partner)) # never yourself
      expect_true(all(groups != groups[partner]))   # never an ingroup
      expect_identical(1:length(groups), sort(partner)) # only one partner per person
    })
  }
}
