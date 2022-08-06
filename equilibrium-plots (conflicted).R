
# Various plots of equilibrium p as function of G and k

# == functions ====
# 
pbar <- function (p, k, G) {
  binoms <- dbinom((k*G):G, G, p)
  # for low probabilities and high G, this gets to 0 and we 
  # then divide by zero
  if (p > 0 && p < 1) binoms <- pmax(binoms, .Machine$double.eps)
  1/G * sum(binoms * (k*G):G) / sum(binoms)
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
  if (rel_fitness(p = 1-1e-4, k = k, G = G, c = c, b = b) > 0) return(1)
  
  rel_fitness_here <- function (p) rel_fitness(p = p, k = k, G = G, c = c, b = b)
  uniroot_result <- uniroot(rel_fitness_here, interval = c(1e-4, 1 - 1e-4))
  uniroot_result$root
}

# equilibrium_p(k = 0.8, G = 10, c = 1, b = 2)

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
  k <- round(10 * k, 1) / 10
  Vectorize(equilibrium_p)(k = k, G = G, c = 1, b = b)
}

df <- expand.grid(k = seq(0.01, 0.99, 0.01), b = seq(1, 10, 0.1), G = c(20, 100, 1000))
df$p <- purrr::pmap_dbl(df, plot_eqm_3d)
#df$p[df$p == 0] <- NA_real_
ggplot(df, aes(x = k, y = b, fill = p)) + 
   geom_raster() + 
   scale_fill_viridis_c() +
   scale_x_continuous(breaks = seq(0.1, 0.9, 0.1)) +
   scale_y_continuous(breaks = 1:10) +
   facet_wrap(vars(G), labeller = label_both) +
   labs(x = "Reciprocation threshold (k)", y = "Benefit/cost (b/c)") +
   theme_minimal() + 
   theme(axis.text.x = element_text(size = 8))

ggsave("heatmap.pdf", bg = "white", width = 10, height = 5)


