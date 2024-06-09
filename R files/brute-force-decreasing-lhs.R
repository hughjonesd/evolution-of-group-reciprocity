
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


# New stuff by Ro'i

q <- function (p, k, G) {
  sum(
    dbinom((k*G):G, G, p)
  )
}

fitnessGR <- function (p, k, G, bc) {
  q(p, k, G) * (pbar(p, k, G) * bc - 1) * q(p, k, G) * (pbar(p, k, G) * q(p, k, G) / p)
}

fitnessS <- function (p, k, G, bc) {
  (q(p, k, G) * (pbar(p, k, G) * bc) * q(p, k, G)) * ((1 - pbar(p, k, G) * q(p, k, G)) / (1 - p))
}

diffFitness <- function (p, k, G, bc) {
  fitnessGR(p, k, G, bc) - fitnessS(p, k, G, bc)
}

# Plots
library(ggplot2)

ps <- seq(0, 100, 1)
psp <- ps / 100

forplot <- data.frame(psp)

plot <- function(k, G) {
  forplot$lhs <- lhs(psp,k,G)
  ggplot(forplot, aes(x = psp, y = lhs)) +
  geom_path(color = "blue", size = 1) +  # Change the line color and size +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +  # Customize x-axis
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +  # Customize y-axis
  geom_hline(yintercept = 0.5, color = "red", linetype = "dashed", size = 1) +  # Add a horizontal line at y = 0.5
  annotate("text", x = Inf, y = 0.5, label = "c/b = 0.5", hjust = 1.1, vjust = -0.5, color = "red", size = 5) +  # Add labe
  theme_minimal() +                      # Use a minimal theme
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center and style the title
    axis.title = element_blank(),               # Style axis titles
    axis.text = element_text(size = 10)                                # Style axis text
  )
}

k_values <- seq(0.5, 0.99, by = 0.1)

for (i in seq_along(k_values)) {
  k <- k_values[i]
  tempplot <- plot(k,10)
  ggsave(filename = paste0("plotk", i, ".png"), plot = tempplot, width = 8, height = 6, dpi = 300, units = "in")
}

G_values <- seq(10, 100, by = 10)

for (i in seq_along(G_values)) {
  G <- G_values[i]
  tempplot <- plot(0.8,G)
  ggsave(filename = paste0("plotG", i, ".png"), plot = tempplot, width = 8, height = 6, dpi = 300, units = "in")
}
