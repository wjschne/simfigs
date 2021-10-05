library(extrafont,quietly = T)
loadfonts("pdf", quiet = T)
loadfonts("win",  quiet = T)
library(lavaan, quietly = T)

# WJSmisc is a custom package that can be installed by 
# the next commented line once
# remotes::install_github("wjschne/WJSmisc")

options(tidyverse.quiet = TRUE)
library(tidyverse)
library(patchwork)
library(mvtnorm)
library(WJSmisc)
library(glue)
library(ggtext)
tar_option_set(packages = c("extrafont", "lavaan", "tidyverse","patchwork", "mvtnorm", "simstandard","ggtext", "glue", "WJSmisc"))

make_model <- function(r_gg = .97, # reliability of g
                       r_ss = .92, # reliability of s
                       r_aa = .92, # reliability of a
                       # Regression coefficients
                       b_a.g = 0.5, # g's effect on a
                       b_a.s = 0.3, # s's effect on a
                       b_s.g = 0.65 # g's effect on s
                       ) {
  r_xx <- c(r_gg, r_ss, r_aa)

# Explained variance
s_by_g <- r_gg * (b_s.g ^ 2)
a_by_g <- r_gg * (b_a.g + b_s.g * b_a.s) ^ 2
a_by_s <- (r_ss - s_by_g) * (b_a.s ^ 2)

# Residual variances
v_s <- r_ss - s_by_g
v_a <- r_aa - a_by_g - a_by_s

# Check estimates
# Asymmetric paths
A <- matrix(c(
  0, 0, 0,
  b_s.g, 0, 0,
  b_a.g, b_a.s, 0
), byrow = T, 3, 3)

iA <- solve(diag(3) - A)

# Symmetric paths
S <- diag(c(r_gg, v_s, v_a))

# Implied covariance of true scores
true_implied_cov <- iA %*% S %*% t(iA)

# Check: diagonal of implied covariance should be equal to reliability coefficients
if (!all(diag(true_implied_cov) == r_xx)) stop("Model not correct.")

# Model for simulation
m <- glue::glue("
s ~ {b_s.g} * g
g ~~ {r_gg} * g
s ~~ {v_s} * s
G ~ 1 * g
G ~~ (1 - {r_gg}) * G
S ~ 1 * s
S ~~ (1 - {r_ss}) * S
a ~ {b_a.g} * g + {b_a.s} * s
a ~~ {v_a} * a
A ~ 1 * a
A ~~ (1 - {r_aa}) * A
")
m
}



make_data <- function(m, n, threshold_g, threshold_s, threshold_a, buffer, meaningful_difference, decision_labels) {
  set.seed(123)
d <- lavaan::simulateData(m, sample.nobs = n) |>
  as_tibble() |>
  select(g, s, a, G, S, A) |>
  mutate(across(.fns = \(x) x * 15 + 100)) |>
  mutate(
    sld_o = (G > threshold_g) &
      (S < threshold_s) &
      (A < threshold_a) &
      (G > S + meaningful_difference) &
      (G > A + meaningful_difference),
    sld_l = (g > threshold_g) &
      (s < threshold_s) &
      (a < threshold_a) &
      (g > s + meaningful_difference) &
      (g > a + meaningful_difference),
    outcome = factor(
      sld_o * 2 + sld_l * 1,
      levels = c(3, 0, 2, 1),
      labels = c("TP", "TN", "FP", "FN")),
    sld_relaxed_l = 
      (g > threshold_g - buffer)  &
      (s < threshold_s + buffer) &
      (a < threshold_a + buffer) &
      (g > s + meaningful_difference - buffer) &
      (g > a + meaningful_difference - buffer), 
    sld_relaxed_o =
      (G > threshold_g - buffer) &
      (S < threshold_s + buffer) &
      (A < threshold_a + buffer) &
      (G > S + meaningful_difference - buffer) &
      (G > S + meaningful_difference - buffer),
    decision_o = factor(sld_o + sld_relaxed_o, 
                      levels = 0:2, 
                      labels = decision_labels),
    decision_l = factor(sld_l + sld_relaxed_l, 
                      levels = 0:2, 
                      labels = decision_labels),
    e_g = G - g,
    e_s = S - s,
    e_a = A - a,
    e_dist = sqrt(e_g ^ 2 + e_s ^ 2 + e_a ^ 2),
    S_dist_max = if_else(S < threshold_s, s - threshold_s, 0),
    A_dist_max = if_else(A < threshold_a, a - threshold_a, 0),
    G_dist_max = if_else(G > threshold_g, threshold_g - g, 0),
    S_dist_min = if_else(s < threshold_s, S - threshold_s, 0),
    A_dist_min = if_else(a < threshold_a, A - threshold_a, 0),
    G_dist_min = if_else(g > threshold_g, threshold_g - G, 0)
         ) 

  d$max_dist <- d %>% select(ends_with("_max")) %>% 
  apply(MARGIN = 1, FUN = max)

  d$min_dist <- d %>% select(ends_with("_min")) %>% 
  apply(MARGIN = 1, FUN = max)

  d
  
}

diagnostic_accuracy <- function(x) {
  n <- unname(sum(x))
    c(
      Sensitivity = unname(x["TP"] / (x["TP"] + x["FN"])),
      Specifivity = unname(x["TN"] / (x["TN"] + x["FP"])),
      PPV = unname(x["TP"] / (x["TP"] + x["FP"])),
      NPV = unname(x["TN"] / (x["TN"] + x["FN"])),
      Prevalence = unname(x["TP"] + x["FN"]) / n,
      `Selection Ratio` = unname(x["TP"] + x["FP"]) / n

      
    )
}

make_max_dist_plot <- function(d, family) {
  d %>%
  dplyr::filter(outcome == "FP" | outcome == "FN") %>%
  mutate(
    max_dist_cut = round(max_dist),
    min_dist_cut = round(min_dist),
    dist = ifelse(outcome == "FP", max_dist_cut, min_dist_cut),
    dist = ifelse(dist < 1, 0, dist),
    dist = ifelse(dist > 9.5, 10, dist)
  ) %>%
  group_by(outcome) %>%
  count(dist) %>%  
  mutate(
    p = n / sum(n),
    cp = cumsum(p),
    cp_label = WJSmisc::prob_label(cp) %>% str_remove("1\\.00")
  ) %>%
  ungroup() %>%
  ggplot(aes(x = dist, y = cp)) +
  geom_step(direction = "mid", color = "royalblue") +
  geom_label(
    aes(label = cp_label),
    data = . %>% dplyr::filter(cp < .995 | cp == 1),
    vjust = -0.2,
    hjust = 0.5,
    color = "gray40",
    family = family,
    label.size = 0,
    label.padding = unit(0.5, "mm"),
    size = ggtext_size(13)
  ) +
  scale_x_continuous(
    "Largest Distance of from the Threshold",
    breaks = seq(-50, 50),
    minor_breaks = NULL,
    labels = \(x) ifelse(x == 0, "<1", ifelse(x == 10, "10+", x))
  ) +
  scale_y_continuous(
    "Cumulative Proportion",
    breaks = seq(0, 1, .1),
    minor_breaks = seq(0, 1, 0.05),
    limits = c(0, 1),
    labels = prob_label
  ) +
  facet_grid(rows = vars(factor(
    outcome,
    levels = c("FP", "FN"),
    labels = c("False\nPositives", "False\nNegatives")
  ))) +
  theme_minimal(base_size = 14, base_family = family) +
  theme(strip.text.y = element_text(angle = 0, color = "gray30"), 
        strip.background = element_rect(color = "gray90", size = 1,
                                        fill = "gray90"),
        panel.background = element_rect(color = "gray90",size = 1),
        panel.grid.major = element_line(size = 0.5)) +
  coord_cartesian(clip = "on", ylim = c(0, 1.05))
  
  ggsave("error_distance_from_threshold.pdf", width = 6.5, height = 6.5, dev = cairo_pdf)
  ggsave("error_distance_from_threshold.png", width = 6.5, height = 6.5, dev = ragg::agg_png)
  
}


make_criteria_plot <- function(d, 
                               threshold_g, 
                               threshold_s, 
                               threshold_a, 
                               buffer, 
                               meaningful_difference, 
                               m_cor, 
                               family, 
                               decision_labels) {
  
  tr_s <- threshold_s + buffer
  ts_s <- threshold_s
  tr_a <- threshold_a + buffer
  ts_a <- threshold_a
  tr_g <- threshold_g - buffer
  ts_g <- threshold_g
  
  # Subsample
d_criteria <- d %>% sample_n(4000)

# Base font size
b_size <- 11

# General and Specific
d_polygon_GS <- bind_rows(
  tibble(
    G = c(160, 160, rep(ts_g, length(seq(40, ts_s))), seq(ts_s, 160)),
    S = c(ts_s, 40, seq(40, ts_s), rep(ts_s, length(seq(ts_s, 160)))),
    decision_o = "SLD-Likely") %>% 
    dplyr::filter(G - S >= meaningful_difference),
  tibble(
    G = c(rep(ts_g, length(seq(40, ts_s))), seq(ts_s, 160)),
    S = c(seq(40, ts_s), rep(ts_s, length(seq(ts_s, 160)))),
    decision_o = "SLD-Possible") %>% 
    dplyr::filter(G - S >= meaningful_difference),
   tibble(
    G = c(seq(160, tr_s), rep(tr_g, length(seq(tr_s, 40)))),
    S = c(rep(tr_s, length(seq(160, tr_s))), seq(tr_s, 40)),
    decision_o = "SLD-Possible") %>% 
    dplyr::filter(G - S >= meaningful_difference - buffer),
  tibble(G = c(rep(tr_g, length(seq(40, tr_s))), seq(tr_g, 160)),
         S = c(seq(40, tr_s), rep(tr_s, length(seq(tr_g, 160)))),
           decision_o = "SLD-Unlikely") %>% 
      dplyr::filter(G - S >= meaningful_difference - buffer) %>% 
    add_row(G = c(160, 40, 40),
            S = c(160, 160, 40),
            decision_o = "SLD-Unlikely")
  )

plot_GS <- ggplot(d_criteria, aes(S, G)) + 
  geom_polygon(data = cor_ellipse(m_cor["G", "S"], 
                                           mean = c(100, 100), 
                                           sd = c(15,15), 
                                           p = .95), 
               aes(x = x, y = y), alpha = .2)  +
  geom_polygon(data = d_polygon_GS,
    aes(fill = decision_o), alpha = 0.3) +
  # geom_point(pch = 16, size = 0.2, alpha = 0.2) +
  geom_text(data = tibble(
    S = 41,
    G = c(92.5, 87.5, 82.5),
    label = rev(decision_labels)),
           size = ggtext_size(b_size),
           hjust = 0,
           family = family, 
           aes(label = label)) +
  scale_x_continuous("Specific Ability",
                     breaks = seq(40, 160, 10),
                     expand = expansion()) +
  scale_y_continuous("General Ability",
                     breaks = seq(40, 160, 10),
                     expand = expansion()) +
  theme_minimal(base_size = b_size, base_family = family) +
  coord_equal(xlim = c(40, 160), ylim = c(40, 160)) +
  scale_fill_viridis_d(begin = .05, end = 0.95) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
plot_GS

# General and Academic

d_polygon_GA <- tibble(
    G = c(
      40,40, 160, 160, 95, 85, 85,
      85, 85, 95, 160, 160, 95, 90, 90,
      90, 90, 95, 160, 160), 
    A = c(40, 160, 160, 90, 90, 80, 40,
          40, 80, 90, 90, 85, 85, 80, 40,
          40, 80, 85, 85, 40),
    decision_o = factor(c(
      rep("SLD-Unlikely", 7),
      rep("SLD-Possible", 8),
      rep("SLD-Likely", 5)),
    levels = decision_labels))

plot_GA <- ggplot(d_criteria, aes(G, A)) + 
  geom_polygon(data = cor_ellipse(m_cor["G", "A"], 
                                           mean = c(100, 100), 
                                           sd = c(15,15), 
                                           p = .95), 
               aes(x = x, y = y), alpha = .2)  +
  geom_polygon(data = d_polygon_GA,
    aes(fill = decision_o), alpha = 0.3) +
  # geom_point(pch = 16, size = 0.2, alpha = 0.2) +
  geom_text(data = tibble(G = c(159, 159, 159),
                          A = c(82.5, 87.5, 92.5),
                          label = rev(decision_labels)),
           size = ggtext_size(b_size),
           hjust = 1,
           family = family, 
           aes(label = label)) +
  scale_x_continuous("General Ability",
                     breaks = seq(40, 160, 10),
                     expand = expansion()) +
  scale_y_continuous("Academic Ability",
                     breaks = seq(40, 160, 10),
                     expand = expansion()) +
  theme_minimal(base_size = b_size, base_family = family) +
  coord_equal(xlim = c(40, 160), ylim = c(40, 160)) +
  scale_fill_viridis_d(begin = .05, end = 0.95) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank())
plot_GA

# Specific and Academic
d_polygon_SA <- tibble(
    S = c(
      40, 40, 160, 160, 90, 90, 40,
      40, 90, 90, 85, 85, 40,
      40, 40, 85, 85), 
    A = c(40, 160, 160, 40, 40, 90, 90,
          90,90,40,40,85,85,
          40,85,85,40),
    decision_o = factor(
      c(rep("SLD-Unlikely", 7),
        rep("SLD-Possible", 6),
        rep("SLD-Likely", 4)), 
      levels = decision_labels))

plot_SA <- ggplot(d_criteria, aes(S, A)) + 
  geom_polygon(data = cor_ellipse(m_cor["S", "A"], 
                                           mean = c(100, 100), 
                                           sd = c(15,15), 
                                           p = .95), 
               aes(x = x, y = y), alpha = .2)  +
  geom_polygon(data = d_polygon_SA,
    aes(fill = decision_o), alpha = 0.3) +
  # geom_point(pch = 16, size = 0.2, alpha = 0.2) +
  geom_text(data = tibble(
    S = 41,
    A = c(82.5, 87.5, 92.5),
    label = rev(decision_labels)),
           size = ggtext_size(b_size),
           hjust = 0,
           family = family, 
           aes(label = label)) +
  scale_x_continuous("Specific Ability",
                     breaks = seq(40, 160, 10),
                     expand = expansion()) +
  scale_y_continuous("Academic Ability",
                     breaks = seq(40, 160, 10),
                     expand = expansion()) +
  theme_minimal(base_size = b_size, base_family = family) +
  coord_equal(xlim = c(40, 160), ylim = c(40, 160)) +
  scale_fill_viridis_d(begin = .05, end = 0.95) +
  theme(legend.position = "none")
plot_SA

d_c <- tibble(Criteria = c("General ability is\naverage or better.",
                    "Specific ability\nis low.",
                    "Academic ability\nis low.",
                    "General ability exceeds\nspecific ability.",
                    "General ability exceeds\nacademic ability."),
       `SLD-Likely` = c("G > 90",
                  "S < 85",
                  "A < 85",
                  "G > S + 10",
                  "G > A + 10"),
       `SLD-Possible` = c("G > 85",
                  "S < 90",
                  "A < 90",
                  "G > S + 5",
                  "G > A + 5")) 

hj <- matrix(c(0, 0.5, 0.5),
             ncol = 3,
             nrow = nrow(d_c),
             byrow = TRUE) %>% as.vector()
x <- matrix(c(0.05, 0.5, .5),
            ncol = 3,
            nrow = nrow(d_c),
            byrow = TRUE) %>% as.vector()

bg_fill_colors <- viridis::viridis(3, begin = .05, end = .95, alpha = .28)[c(1,3,2)]

bg_fill <- matrix(bg_fill_colors, 
                  ncol = 3, 
                  nrow = nrow(d_c),
                  byrow = TRUE) %>% as.vector()

table_theme <-
  gridExtra::ttheme_default(
    base_size = 7.71,
    base_family = family,
    core = list(
      fg_params = list(hjust = hj, x = x),
      bg_params = list(fill = bg_fill),
      padding = grid::unit.c(unit(8.3, "mm"), unit(7.15, "mm"))
    ),
    colhead = list( 
      fg_params = list(
        hjust = c(0, 0.5, 0.5),
        x = c(0.05, 0.5, .5)
      ),
      bg_params = list(fill = bg_fill_colors)
    )
  )


width <- 2.82
plot_GS  +
  gridExtra::tableGrob(d_c,theme = table_theme, rows = NULL) +
  plot_SA +
  plot_GA +
  plot_layout(ncol = 2,
              nrow = 2,
              heights = unit(c(width,width), "in"),
              widths = unit(c(width,width), "in"))


ggsave("criteria.png", width = 6.5, height = 6.5, device = ragg::agg_png)
ggsave("criteria.pdf", width = 6.5, height = 6.5, device = cairo_pdf)
  
}

make_ppvnpv_plot <- function(rxx = 0.92, threshold = 85) {

  x <- seq(40, 160, 0.1)




  crossing(x = x, rxx = rxx, threshold = threshold) %>%
  mutate(see = 15 * sqrt(rxx - rxx ^ 2),
         mu = (x - 100) * rxx + 100,
         p = pnorm(threshold, mu, see),
         d = dnorm(x, 100, 15),
         dp = d * p) %>%
  mutate(rxx_label = paste0("italic(r[xx])=='", WJSmisc::prob_label(rxx), "'")) %>%
  ggplot(aes(x, d)) +
  geom_ribbon(color = NA, fill = "#5DC863", aes(ymin = dp, ymax = d, alpha = x <= threshold)) +
  geom_area(aes(y = dp, alpha = x < threshold), fill = "#3B528B", color = NA) +
  theme_classic(base_family = "Roboto Condensed", base_size = 20) +
  geom_vline(aes(xintercept = threshold), color = "white") +
  # facet_grid(rows = vars(rxx_label), labeller = label_parsed) +
  scale_x_continuous("Observed Score", breaks = seq(40, 160, 5),
                     minor_breaks = seq(40, 160, 5),
                     expand = expansion()) +
  scale_y_continuous(NULL, expand = expansion(mult = c(0,.093))) +
  scale_alpha_manual(values = c(0.5, .8)) +
  ggeasy::easy_remove_y_axis() +
  ggeasy::easy_remove_legend() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(3,5,2,2), "mm"),
        axis.line.x = element_line(size = 0.5, color = "gray20"),
        axis.text.x = element_text(color = c("gray20", NA, NA)),
        axis.ticks = element_line(size = c(.50, 0.25, 0.25), color = "gray20"),
        strip.text.y = element_text(angle = 0)) +
  geom_text(data = tibble(
    x = c(86, 84, 84, 86),
    y = c(0.011, 0.011, 0.001, 0.001),
    l = c("True Negatives", "False Positives", "True Positives", "False Negatives"),
    hjust = c(0, 1, 1, 0)),
    color = "gray20",
    aes(y = y, label = l, hjust = hjust),
    family = "Roboto Condensed",
    size = WJSmisc::ggtext_size(20)) +
    geom_tile(data = tibble(x = c(59, 65, 71, 129, 135, 141),
                            y = c(.0208, .025, .0208, .0208, .025, .0208),
                            fill = c("#5DC863", "#5DC863", "#3B528B", "#3B528B", "#3B528B", "#5DC863")),
              aes(x = x, y = y, fill = fill, alpha = x < threshold), width = 8, height = .003) + 
    scale_fill_identity() + 
    scale_color_identity() + 
    annotate(geom = "segment", x = c(54, 124), xend = c(76, 146), y = (.0208 + .025) / 2, yend = (.0208 + .025) / 2, size = .25, color = "gray30") + 
    geom_text(data = tibble(x = c(59, 65, 71, 129, 135, 141),
                            y = c(.0208, .025, .0208, .0208, .025, .0208),
                            l = c(.13, .13, .03, .82, .82, .02)),
              aes(x = x, y = y, label = WJSmisc::prob_label(l)),
              size = WJSmisc::ggtext_size(18),
              family = "Roboto Condensed",
              color = "gray10"
              ) + 
    geom_text(data = tibble(x = c(65, 135)), aes(x =x, y = .0208, label = "+"), 
             size = WJSmisc::ggtext_size(18), 
             family = "Roboto Condensed",
             color = "gray30") + 
    geom_text(data = tibble(x = c(51.5, 78.5, 121.5, 148.5, 46, 83.5, 116, 153.5), 
                            y = (.0208 + .025) / 2, 
                            label = c(rep("=", 4), "PPV", ".80", "NPV", ".97")),
              aes(x = x, y = y, label = label),
             size = WJSmisc::ggtext_size(20), 
             family = "Roboto Condensed",
             color = "gray30") + 
    geom_text(data = tibble(x = c(59, 65, 71, 129, 135, 141),
                            y = c(.0208, .025, .0208, .0208, .025, .0208),
                            color = c("#5DC863", "#5DC863", "#3B528B", "#3B528B", "#3B528B", "#5DC863"),
                            l = c("True\nPositives", "True\nPositives", "False\nPositives",
                                  "True\nNegatives", "True\nNegatives", "False\nNegatives"),
                            vjust = c(2.06, -1.08, 2.06, 2.06, -1.08,2.06)),
              aes(x = x, y = y, color = color, alpha = x < threshold, 
                  label = l, vjust = vjust), lineheight = .9,
              size = WJSmisc::ggtext_size(13), 
             family = "Roboto Condensed") + 
    annotate("text", x = 126, y = 0.011, 
             label = paste0("Reliability = ",WJSmisc::prob_label(rxx),"\nThreshold = ", threshold), 
             size = ggtext_size(20), hjust = 0, family = "Roboto Condensed", lineheight = .95, color = "gray30")


ggsave("ppvnpv.png", height = 4.9, width = 6.5, device = ragg::agg_png)
ggsave("ppvnpv.pdf", height = 4.9, width = 6.5, device = cairo_pdf)

}


options(scipen = 999)


conditional_ppv <- function(General = 100,
                            Specific = 85,
                            Academic = 70,
                            # Reliability
                            r_gg = 0.97,
                            r_ss = 0.92,
                            r_aa = 0.96,
                            # Thresholds
                            threshold_s = 85,
                            threshold_a = 85,
                            threshold_g = 90,
                            buffer = 5,
                            meaningful_difference = 10,
                            # Regression effects
                            b_a.g = 0.5, # g's effect on a
                            b_a.s = 0.3, # s's effect on a
                            b_s.g = 0.65, # g's effect on s
                            myfont = "Roboto Condensed",
                            filename,
                            d) {
  
  # Subsample
  d_criteria <- d %>% sample_n(2000) 
  
  # Base font size
  b_size <- 11
  
  # Observed Scores
  x_SS <- c(G = General, S = Specific, A = Academic)
  
  # Variable Names
  v_true <- c("g", "s", "a")
  v_observed <- c("G", "S", "A")
  v_names <- c(v_true, v_observed)
  decision_labels <- c("SLD-Unlikely",
                       "SLD-Possible",
                       "SLD-Likely")
  
  # Reliability
  r_xx <- c(r_gg, r_ss, r_aa)
  
  
  # Explained variance
  s_by_g <- r_gg * (b_s.g ^ 2)
  a_by_g <- r_gg * (b_a.g + b_s.g * b_a.s) ^ 2
  a_by_s <- (r_ss - s_by_g) * (b_a.s ^ 2)
  
  # Residual variances
  v_s <- r_ss - s_by_g
  v_a <- r_aa - a_by_g - a_by_s
  
  
  # Model for simulation
  m <- glue::glue(
    "
s ~ {b_s.g} * g
g ~~ {r_gg} * g
s ~~ {v_s} * s
G ~ 1 * g
G ~~ (1 - {r_gg}) * G
S ~ 1 * s
S ~~ (1 - {r_ss}) * S
a ~ {b_a.g} * g + {b_a.s} * s
a ~~ {v_a} * a
A ~ 1 * a
A ~~ (1 - {r_aa}) * A
"
  )
  
  # Model-implied Covariance Matrix
  m_cov <- fitted(sem(m))$cov[v_names, v_names]
  m_cor <- cov2cor(m_cov)
  cov_all <- m_cov * 225
  # True score covariance
  cov_true <- cov_all[v_true, v_true]
  # Observed score covariance
  cov_observed <- cov_all[v_observed, v_observed]
  # Cross covariances
  cov_true_observed <-  cov_all[v_true, v_observed]
  
  # Conditional means
  mu_conditional <-
    100 + 
    cov_true_observed %*% 
    solve(cov_observed) %*% 
    (x_SS - 100) %>% 
    as.numeric()
  
  # Conditional covariance
  cov_conditional <-
    cov_true - 
    cov_true_observed %*% 
    solve(cov_observed) %*% 
    t(cov_true_observed)
  
  # Difference weights
  w_difference <- (
    "
variable	gs   	ga   	GS   	GA
g       	1.00 	1.00 	0.00 	0.00
s       	-1.00	0.00 	0.00 	0.00
a       	0.00 	-1.00	0.00 	0.00
G       	0.00 	0.00 	1.00 	1.00
S       	0.00 	0.00 	-1.00	0.00
A       	0.00 	0.00 	0.00 	-1.00"
  ) %>%
    I() %>% 
    readr::read_tsv(file = ., show_col_types = F) %>%
    tibble::column_to_rownames("variable") %>%
    select(-GS,-GA) %>%
    as.matrix()
  
  # Full weight matrix
  w <- cbind(diag(6), w_difference) %>%
    `colnames<-`(c(rownames(w_difference),
                   colnames(w_difference)))
  
  # Conditional Correlations
  cor_conditional <- cov2cor(cov_conditional)
  
  # Covariance matrix of all scores and true difference scores
  big_sigma <- composite_covariance(cov_all, w)
  # Conditional covariance of all true scores and difference scores
  cond_mu_sigma <-
    condMVNorm::condMVN(
      mean = c(rep(100, 6), rep(0, 2)) %>% `names<-`(rownames(big_sigma)),
      sigma = big_sigma,
      dependent.ind = c("g", "s", "a", "gs", "ga"),
      given.ind = toupper(c("g", "s", "a")),
      X.given = x_SS,
      check.sigma = F
    )
  
  # Percent meeting strict criteria
  p_strict <- pmvnorm(
    lower = c(
      threshold_g,-Inf,-Inf,
      meaningful_difference,
      meaningful_difference
    ),
    upper = c(Inf,
              threshold_s,
              threshold_a,
              Inf,
              Inf),
    mean = cond_mu_sigma$condMean,
    sigma = cond_mu_sigma$condVar,
    keepAttr = F
  )
  
  # Percent meeting relaxed criteria
  p_relaxed <- pmvnorm(
    lower = c(
      threshold_g - buffer,-Inf,-Inf,
      meaningful_difference - buffer,
      meaningful_difference - buffer
    ),
    upper = c(Inf,
              threshold_s + buffer,
              threshold_a + buffer,
              Inf,
              Inf),
    mean = cond_mu_sigma$condMean,
    sigma = cond_mu_sigma$condVar,
    keepAttr = F
  )

  
  
  
  # Simulated true scores conditioned on observed scores
  d_sa <-
    rmvnorm(n = 100, 
            mean = mu_conditional %>% 
              `names<-`(c("g", "s", "a")), 
            sigma = cov_conditional) %>%
    as_tibble() %>%
    mutate(
      sld_l = (g > threshold_g) &
        (s < threshold_s) &
        (a < threshold_a) &
        (g > s + meaningful_difference) &
        (g > a + meaningful_difference),
      sld_relaxed_l =
        (g > threshold_g - buffer)  &
        (s < threshold_s + buffer) &
        (a < threshold_a + buffer) &
        (g > s + meaningful_difference - buffer) &
        (g > a + meaningful_difference - buffer),
      decision_l = factor(
        sld_l + sld_relaxed_l,
        levels = 0:2,
        labels = decision_labels
      )
    )

  tr_s <- threshold_s + buffer
  ts_s <- threshold_s
  tr_a <- threshold_a + buffer
  ts_a <- threshold_a
  tr_g <- threshold_g - buffer
  ts_g <- threshold_g
  # General and Specific ----
# d_polygon_GS <- tibble(
#     G = c(
#       40,40, 160, 160, ta_s, ta_g, ta_g,
#       ta_g, ta_g, ta_s, 160, 160, ta_s, tr_g, tr_g,
#       tr_g, tr_g, ta_s, 160, 160, ta_s, ts_g, ts_g,
#       ts_g, ts_g, ta_s, 160, 160), 
#     S = c(40, 160, 160, ta_s, ta_s, ta_g, 40,
#           40, ta_g, ta_s, ta_s, tr_s, tr_s, ta_g, 40,
#           40, ta_g, tr_s, tr_s, ts_s, ts_s, ta_g, 40,
#           40, ta_g, ts_s, ts_s, 40),
#     decision_o = factor(
#       c(rep("Not SLD", 7),
#         rep("SLD-Adjacent", 8),
#         rep("SLD-Relaxed", 8),
#         rep("SLD-Strict", 5)), 
#       levels = decision_labels))
  
d_polygon_GS <- bind_rows(
  tibble(
    G = c(160, 160, rep(ts_g, length(seq(40, ts_s))), seq(ts_s, 160)),
    S = c(ts_s, 40, seq(40, ts_s), rep(ts_s, length(seq(ts_s, 160)))),
    decision_o = "SLD-Likely") %>% 
    dplyr::filter(G - S >= meaningful_difference),
  tibble(
    G = c(rep(ts_g, length(seq(40, ts_s))), seq(ts_s, 160)),
    S = c(seq(40, ts_s), rep(ts_s, length(seq(ts_s, 160)))),
    decision_o = "SLD-Possible") %>% 
    dplyr::filter(G - S >= meaningful_difference),
   tibble(
    G = c(seq(160, tr_s), rep(tr_g, length(seq(tr_s, 40)))),
    S = c(rep(tr_s, length(seq(160, tr_s))), seq(tr_s, 40)),
    decision_o = "SLD-Possible") %>% 
    dplyr::filter(G - S >= meaningful_difference - buffer),
  tibble(G = c(rep(tr_g, length(seq(40, tr_s))), seq(tr_g, 160)),
         S = c(seq(40, tr_s), rep(tr_s, length(seq(tr_g, 160)))),
           decision_o = "SLD-Unlikely") %>% 
      dplyr::filter(G - S >= meaningful_difference - buffer) %>% 
    add_row(G = c(160, 40, 40),
            S = c(160, 160, 40),
            decision_o = "SLD-Unlikely")
  )


plot_GS <- ggplot(d_criteria, aes(S, G)) + 
  geom_polygon(data = cor_ellipse(m_cor["G", "S"],
                                  mean = c(100, 100),
                                  sd = c(15,15),
                                  p = .95),
               aes(x = x, y = y), alpha = .2)  +
  geom_polygon(data = d_polygon_GS,
    aes(fill = decision_o), alpha = 0.3) +
  # geom_point(pch = 16, size = 0.2, alpha = 0.2) +
  geom_text(data = tibble(
    S = 41,
    G = c(tr_g - buffer / 2, ts_g - buffer / 2, ts_g + buffer / 2),
    label = (decision_labels)),
           size = ggtext_size(b_size),
           hjust = 0,
           family = myfont, 
           aes(label = label)) +
  scale_x_continuous("Specific Ability",
                     breaks = seq(40, 160, 10),
                     expand = expansion()) +
  scale_y_continuous("General Ability",
                     breaks = seq(40, 160, 10),
                     expand = expansion()) +
  theme_minimal(base_size = b_size, base_family = myfont) +
  coord_equal(xlim = c(40, 160), ylim = c(40, 160)) +
  scale_fill_viridis_d(begin = .05, end = 0.95) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
# plot_GS

# General and Academic ----

d_polygon_GA <- bind_rows(
  tibble(
    G = c(160, 160, rep(ts_g, length(seq(40, ts_a))), seq(ts_a, 160)),
    A = c(ts_a, 40, seq(40, ts_a), rep(ts_a, length(seq(ts_a, 160)))),
    decision_o = "SLD-Likely") %>% 
    dplyr::filter(G - A >= meaningful_difference),
  tibble(
    G = c(rep(ts_g, length(seq(40, ts_a))), seq(ts_a, 160)),
    A = c(seq(40, ts_a), rep(ts_a, length(seq(ts_a, 160)))),
    decision_o = "SLD-Possible") %>% 
    dplyr::filter(G - A >= meaningful_difference),
   tibble(
    G = c(seq(160, tr_a), rep(tr_g, length(seq(tr_a, 40)))),
    A = c(rep(tr_a, length(seq(160, tr_a))), seq(tr_a, 40)),
    decision_o = "SLD-Possible") %>% 
    dplyr::filter(G - A >= meaningful_difference - buffer),
  tibble(G = c(rep(tr_g, length(seq(40, tr_a))), seq(tr_g, 160)),
         A = c(seq(40, tr_a), rep(tr_a, length(seq(tr_g, 160)))),
           decision_o = "SLD-Unlikely") %>% 
      dplyr::filter(G - A >= meaningful_difference - buffer) %>% 
    add_row(G = c(160, 40, 40),
            A = c(160, 160, 40),
            decision_o = "SLD-Unlikely")
  )



plot_GA <- ggplot(d_criteria, aes(G, A)) + 
  geom_polygon(data = cor_ellipse(m_cor["G", "A"], 
                                           mean = c(100, 100), 
                                           sd = c(15,15), 
                                           p = .95), 
               aes(x = x, y = y), alpha = .2)  +
  geom_polygon(data = d_polygon_GA,
    aes(fill = decision_o), alpha = 0.3) +
  # geom_point(pch = 16, size = 0.2, alpha = 0.2) +
  geom_text(data = tibble(G = c(159, 159, 159),
                          A = c(tr_a + buffer / 2, 
                                ts_a + buffer / 2, 
                                ts_a - buffer / 2),
                          label = (decision_labels)),
           size = ggtext_size(b_size),
           hjust = 1,
           family = myfont, 
           aes(label = label)) +
  scale_x_continuous("General Ability",
                     breaks = seq(40, 160, 10),
                     expand = expansion()) +
  scale_y_continuous("Academic Ability",
                     breaks = seq(40, 160, 10),
                     expand = expansion()) +
  theme_minimal(base_size = b_size, base_family = myfont) +
  coord_equal(xlim = c(40, 160), ylim = c(40, 160)) +
  scale_fill_viridis_d(begin = .05, end = 0.95) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank())




# Specific and Academic ----
d_polygon_SA <- tibble(
    S = c(
      40, 40, 160, 160, tr_s, tr_s, 40,
      40, tr_s, tr_s, ts_s, ts_s, 40,
      40, 40, ts_s, ts_s), 
    A = c(40, 160, 160, 40, 40, tr_a, tr_a,
          tr_a,tr_a,40,40,ts_a,ts_a,
          40,ts_a,ts_a,40),
    decision_o = factor(
      c(rep("SLD-Unlikely", 7),
        rep("SLD-Possible", 6),
        rep("SLD-Likely", 4)), 
      levels = decision_labels))



plot_SA <- ggplot(d_criteria, aes(S, A)) + 
  geom_polygon(data = cor_ellipse(m_cor["S", "A"], 
                                           mean = c(100, 100), 
                                           sd = c(15,15), 
                                           p = .95), 
               aes(x = x, y = y), alpha = .2)  +
  geom_polygon(data = d_polygon_SA,
    aes(fill = decision_o), alpha = 0.3) +
  # geom_point(pch = 16, size = 0.2, alpha = 0.2) +
  geom_text(data = tibble(
    S = 41,
    A = c(tr_a + buffer / 2, 
                                ts_a + buffer / 2, 
                                ts_a - buffer / 2),
    label = decision_labels),
           size = ggtext_size(b_size),
           hjust = 0,
           family = myfont, 
           aes(label = label)) +
  scale_x_continuous("Specific Ability",
                     breaks = seq(40, 160, 10),
                     expand = expansion()) +
  scale_y_continuous("Academic Ability",
                     breaks = seq(40, 160, 10),
                     expand = expansion()) +
  theme_minimal(base_size = b_size, base_family = myfont) +
  coord_equal(xlim = c(40, 160), ylim = c(40, 160)) +
  scale_fill_viridis_d(begin = .05, end = 0.95) +
  theme(legend.position = "none")

  case_SA <-
    plot_SA +
    geom_point(
      data = d_sa,
      aes(x = s, y = a),
      alpha = .3,
      pch = 16,
      size = .1
    )  +
    scale_color_viridis_d() +
    geom_polygon(
      d = cor_ellipse(
        r = cor_conditional[2, 3],
        mean = mu_conditional[c(2, 3)],
        sd = sqrt(diag(cov_conditional)[c(2, 3)]),
        p = .95
      ),
      aes(x = x, y = y),
      fill = "firebrick4",
      alpha = .2
    ) +
    annotate("point", x = x_SS["S"], y = x_SS["A"]) +
    annotate(
      "text",
      hjust = 1.1,
      vjust = 1.1,
      family = myfont,
      x = x_SS["S"],
      y = x_SS["A"],
      label = paste0("(", x_SS["S"], ",", x_SS["A"], ")")
    )
  case_GA <- plot_GA +
    geom_point(
      data = d_sa,
      aes(x = g, y = a),
      alpha = .3,
      pch = 16,
      size = .1
    )   +
    scale_color_viridis_d() +
    geom_polygon(
      d = cor_ellipse(
        r = cor_conditional[1, 3],
        mean = mu_conditional[c(1, 3)],
        sd = sqrt(diag(cov_conditional)[c(1, 3)]),
        p = .95
      ),
      aes(x = x, y = y),
      fill = "firebrick4",
      alpha = .2
    ) +
    annotate("point", x = x_SS["G"], y = x_SS["A"]) +
    annotate(
      "text",
      hjust = -.1,
      vjust = 1.1,
      family = myfont,
      x = x_SS["G"],
      y = x_SS["A"],
      label = paste0("(", x_SS["G"], ",", x_SS["A"], ")")
    )
  
  
  
  case_GS <-
    plot_GS +
    geom_point(
      data = d_sa,
      aes(x = s, y = g),
      alpha = .3,
      pch = 16,
      size = 0.1
    ) +
    scale_color_viridis_d() +
    geom_polygon(
      d = cor_ellipse(
        r = cor_conditional[2, 1],
        mean = mu_conditional[c(2, 1)],
        sd = sqrt(diag(cov_conditional)[c(2, 1)]),
        p = .95
      ),
      aes(x = x, y = y),
      fill = "firebrick4",
      alpha = .2
    ) +
    annotate(
      "text",
      hjust = 1.1,
      vjust = -.1,
      family = myfont,
      x = x_SS["S"],
      y = x_SS["G"],
      label = paste0("(", x_SS["S"], ",", x_SS["G"], ")")
    ) +
    annotate("point", x = x_SS["S"], y = x_SS["G"])
  
  grob_outcome <- tibble(
    `Conditional\nOutcome` = rev(decision_labels),
    `Outcome\nProportion` = c(
      p_strict,
      p_relaxed - p_strict,
      1 - p_relaxed
    ),
    `Cumulative\nProportion` = cumsum(`Outcome\nProportion`)
  ) %>%
    mutate(across(
      where(is.numeric),
      \(x) WJSmisc::prob_label(x, digits = 2, max_digits = 6)
    )) %>%
    gridExtra::tableGrob(
      rows = NULL,
      theme = gridExtra::ttheme_default(
        base_size = 11.5,
        base_family = myfont,
        core = list(
          fg_params = list(hjust = c(rep(0, 3), rep(0.5, 6)),
                           x = c(rep(0.07, 3), rep(0.5, 6))),
          bg_params = list(fill = rep(
            c("#F6F7BA", "#C1E9D8", "#C1D3DE"), 3
          )),
          padding = grid::unit.c(unit(5.0, "mm"), unit(16.2, "mm"))
        ),
        colhead = list(
          fg_params = list(
            hjust = c(0, 0.5, 0.5),
            x = c(0.07, 0.5, 0.5)
          ),
          bg_params = list(fill = c("white"))
        )
      )
    )
  
  
  case_GS + grob_outcome + case_SA + case_GA
  
  ggsave(filename = paste0(filename,".pdf"), width = 6.5, height = 6.5, device = cairo_pdf)
  ggsave(filename = paste0(filename,".png"), width = 6.5, height = 6.5, device = ragg::agg_png)
}



