
# test whether (pbar - p)/(1-p) is decreasing for many
# values of G and k = R/G

pbar <- function (p, k, G) {
  1/G * sum(
          dbinom((k*G):G, G, p) * (k*G):G
        ) /
        sum(
          dbinom((k*G):G, G, p)
        )
}
pbar <- Vectorize(pbar)

lhs <- function (p, k, G) {
  (pbar(p, k, G) - p) / (1 - p)
}


max_G <- 1000
ps <- seq(0, 1, 100)
pb <- progress::progress_bar$new(total = max_G * (max_G + 1)/2, clear = FALSE)

for (G in 1:max_G) for (R in 1:G) {
  pb$tick()
  k <- R/G
  lhs_vals <- lhs(ps, k = k, G = G)
  diff_lhs <- diff(lhs_vals)
  lhs_increasing <- diff_lhs >= 0
  if (any(lhs_increasing)) {
    starts <- lhs_vals[lhs_increasing]
    ends   <- lhs_vals[lhs_increasing + 1]
    p_starts <- ps[lhs_increasing]
    p_ends   <- ps[lhs_increasing + 1]
    stop("LHS of (2) was not decreasing between these values:\n", 
           sprintf("p = %s (LHS = %s) and p = %s (LHS = %s)", p_starts, starts, 
                     p_ends, ends)
         )
  }
}
