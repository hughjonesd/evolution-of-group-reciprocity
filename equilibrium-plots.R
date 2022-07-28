
# Various plots of equilibrium p as function of G and k

# == functions ====
# 
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

rel_fitness <- function (p, k, G, c, b) {
  lhs(p, k, G) - c/b
}

equilibrium_p <- function (k, G, c, b) {
  if (rel_fitness(p = 1e-4, k = k, G = G, c = c, b = b) < 0) return(0)
  
  rel_fitness_here <- function (p) rel_fitness(p = p, k = k, G = G, c = c, b = b)
  uniroot_result <- uniroot(rel_fitness_here, interval = c(1e-4, 1 - 1e-4))
  uniroot_result$root
}

equilibrium_p(k = 0.8, G = 10, c = 1, b = 2)

# == ggplot version ====

library(ggplot2)
eqm_function  <- function (k) {
  f <- function (x) Vectorize(equilibrium_p)(k = k, G = 10, c = 1, b = x)
  
  f
}

line_width <- 1.25
ggplot() + xlim(1, 10) + 
  geom_function(
                fun = eqm_function(k = 0.4), 
                aes(color = "k = 0.4"),
                size = line_width
                ) + 
  geom_function(
                fun = eqm_function(k = 0.6), 
                aes(color = "k = 0.6"),
                size = line_width,
                linetype = 2
                ) +
  geom_function(
                fun = eqm_function(k = 0.8), 
                aes(color = "k = 0.8"),
                size = line_width,
                linetype = 4
                ) +
  scale_x_continuous(breaks = 1:10, limits = c(1, 10), ) + 
  theme_minimal() + theme(panel.grid.minor = element_blank()) + 
  labs(color = "", x = "Benefit/cost ratio (b/c)", 
       y = "Equilibrium proportion of group reciprocators (p)",
       main = "Group reciprocity")
 

# == 3d version ====
library(rgl)
plot_eqm_3d <- function (k, b) {
  k <- round(2 * k, 1) /2
  Vectorize(equilibrium_p)(k = k, G = 20, c = 1, b = b)
}

rgl::plot3d(plot_eqm_3d, xlim = c(0.05, 0.95), ylim = c(2, 10), 
            zlab = "Equilibrium p", col = hcl.colors)
rgl::snapshot3d("eqm3d.png")

# == heatmap version ====

library(ggplot2)
plot_eqm_3d <- function (k, b, G) {
  k <- round(2 * k, 1) /2
  Vectorize(equilibrium_p)(k = k, G = 20, c = 1, b = b)
}

for (G in c(20, 100, 200)) {
  df <- expand.grid(k = seq(0.05, 0.95, 0.05), b = 1:10)
  df$p <- purrr::map2_dbl(df$k, df$b, plot_eqm_3d, G = G)
  ggp <- ggplot(df, aes(x = k, y = b, fill = p)) + geom_tile() + theme_minimal() + 
           scale_fill_viridis_c() +
           scale_x_continuous(breaks = seq(0.1, 0.9, 0.1)) +
           scale_y_continuous(breaks = 1:10) +
           labs(x = "Reciprocation threshold (k)", y = "Benefit/cost (b/c)",
                  title = "Equilibrium proportion of reciprocators", 
                  subtitle = sprintf("G = %s", G))
  ggsave(sprintf("heatmap-G-%s.jpeg", G), plot = ggp, bg = "white")
}
