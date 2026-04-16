library(ggplot2)
library(scales)
library(gridExtra)
library(grid)
library(dplyr)
library(tidyr)
#Run after Working_Test.R
source("~/Documents/DISSERTATION/Project/Data and Utils.R")
source("~/Documents/DISSERTATION/Project/Deterministic.R")
source("~/Documents/DISSERTATION/Project/Shock.R")
source("~/Documents/DISSERTATION/Project/Metrics.R")
source("~/Documents/DISSERTATION/Project/Scenario.R")
source("~/Documents/DISSERTATION/Project/Working Test.R")




#Deterministic fitted outputs
clm_det <- clm(X_upper)
dcl_det_par <- dcl.estimation(X_upper, N_upper, Tables = FALSE, adj = 1)
dcl_det_fit <- dcl.predict(dcl_det_par, N_upper, Model = 2, Tail = FALSE, Tables = FALSE)
bdcl_det_par <- bdcl.estimation(X_upper, N_upper, I_upper, Tables = FALSE, adj = 1)
bdcl_det_fit <- dcl.predict(bdcl_det_par, N_upper, Model = 2, Tail = FALSE, Tables = FALSE)

make_cy_df <- function(truth, clm_fit, dcl_fit, bdcl_fit) {
  t_cy <- outstanding_by_cy(truth)
  data.frame(
    CY = as.numeric(names(t_cy)),
    Truth = as.numeric(t_cy),
    CLM = as.numeric(outstanding_by_cy(clm_fit)),
    DCL = as.numeric(outstanding_by_cy(dcl_fit)),
    BDCL = as.numeric(outstanding_by_cy(bdcl_fit))
  )
}

cy_det_df <- make_cy_df(X_det, clm_det$triangle.hat,
                        dcl_det_fit$Xtotal, bdcl_det_fit$Xtotal)
cy_sk_df  <- make_cy_df(X_sk,  clm_sk$triangle.hat,
                        dcl_test_fit$Xtotal, bdcl_test_fit$Xtotal)


th <- theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA))

cols <- c(CLM = "brown2", DCL = "skyblue3", BDCL = "springgreen3") 
sub <- bquote(y == .(test_y) ~ "," ~ omega == .(test_omega))
cys <- cy_sk_df$CY





################################################################################

#Reserve Bias reference scenario, Figure 2
total_df <- data.frame(
  model = c("Truth", "CLM", "DCL", "BDCL"),
  total = c(
    sum(outstanding_by_cy(X_sk)),
    sum(outstanding_by_cy(clm_sk$triangle.hat)),
    sum(outstanding_by_cy(dcl_test_fit$Xtotal)),
    sum(outstanding_by_cy(bdcl_test_fit$Xtotal))
  )
) %>%
  mutate(
    model = factor(model, levels = c("Truth","CLM","DCL","BDCL")),
    bias = (total - total[model == "Truth"]) / total[model == "Truth"],
    bar_fill = case_when(
      model == "Truth" ~ "black",
      bias > 0 ~ "steelblue",  #over-reserve
      bias < 0 ~ "brown2"   #under-reserve
    )
  )

truth_val <- total_df$total[total_df$model == "Truth"]

rbias <- ggplot(total_df, aes(model, total, fill = bar_fill)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = truth_val, linetype = "dashed", colour = "black") +
  geom_text(aes(label = paste0(ifelse(bias == 0, "",
                                      ifelse(bias > 0, "+", "")),
                               round(bias * 100, 1), "%")),
            vjust = -0.5, size = 4) +
  scale_fill_identity() +
  scale_y_continuous(labels = label_comma(scale = 1e-6, suffix = "M"),
                     expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Total Reserve by Model vs True Outstanding Reserves",
    subtitle = bquote(y == .(test_y) ~ "," ~ omega == .(test_omega)),
    x = NULL,
    y = "Total Outstanding Reserve"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position  = "none")
rbias

################################################################################

#Calendar-year relative error reference scenario, Figure 3
err_df <- rbind(
  data.frame(CY = cys, err = as.numeric(cy_rel_err(X_sk, clm_sk$triangle.hat)), 
             model = "CLM"),
  data.frame(CY = cys, err = as.numeric(cy_rel_err(X_sk, dcl_test_fit$Xtotal)), 
             model = "DCL"),
  data.frame(CY = cys, err = as.numeric(cy_rel_err(X_sk, bdcl_test_fit$Xtotal)), 
             model = "BDCL")
)


cyre <- ggplot(err_df, aes(CY, err, colour = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line() + geom_point(size = 2) +
  scale_colour_manual(values = cols) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Calendar-Year Relative Error", subtitle = sub,
       x = "Calendar Year", y = "(pred - truth) / truth") + th
cyre

################################################################################

#Wasserstein Distance Reference Scenario, Figure 4
make_pmf <- function(tri, label) {
  cy_vals <- outstanding_by_cy(tri)
  total   <- sum(cy_vals, na.rm = TRUE)
  data.frame(
    CY = as.numeric(names(cy_vals)),
    prob = as.numeric(cy_vals) / total,
    model = label
  )
}

pmf_truth <- make_pmf(X_sk, "Truth")
pmf_clm <- make_pmf(clm_sk$triangle.hat, "CLM")
pmf_dcl <- make_pmf(dcl_test_fit$Xtotal, "DCL")
pmf_bdcl <- make_pmf(bdcl_test_fit$Xtotal, "BDCL")

truth_cdf <- pmf_truth %>%
  mutate(cdf_truth = cumsum(prob)) %>%
  select(CY, cdf_truth)

gap_df <- bind_rows(
  pmf_clm %>% mutate(model = "CLM"),
  pmf_dcl %>% mutate(model = "DCL"),
  pmf_bdcl %>% mutate(model = "BDCL")
) %>%
  group_by(model) %>%
  mutate(cdf_model = cumsum(prob)) %>%
  ungroup() %>%
  left_join(truth_cdf, by = "CY") %>%
  mutate(gap = abs(cdf_model - cdf_truth),
         model = factor(model, levels = c("CLM","DCL","BDCL")))

emd_labels <- gap_df %>%
  group_by(model) %>%
  summarise(emd_val = sum(gap), .groups = "drop") %>%
  mutate(label = paste0("Wasserstein Distance = ",
                        formatC(emd_val, format = "f", digits = 4)))

y_max <- max(gap_df$gap) * 1.15

emd_ref <- ggplot(gap_df, aes(CY, gap)) +
  geom_col(aes(fill = model), width = 0.8, alpha = 0.7) +
  geom_text(data = emd_labels,
            aes(x = 28, y = y_max * 0.9, label = label),
            size = 4, fontface = "bold", inherit.aes = FALSE) +
  facet_wrap(~model, nrow = 1) +
  scale_fill_manual(values = c(CLM = "brown2", DCL = "skyblue3",
                               BDCL = "springgreen3"), guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 0.1),
                     limits = c(0, y_max)) +
  labs(title = "1-Wasserstein Distance Decomposition: |Model CDF - Truth CDF| by Calendar Year",
       subtitle = bquote("Sum of bar heights = Wasserstein Distance; " ~
                           y == .(test_y) ~ "," ~ omega == .(test_omega)),
       x = "Calendar-Year", y = "|CDF difference|") +
  th +
  theme(strip.text = element_text(face = "bold"))
emd_ref


################################################################################

#gamma_i reference scenario, Figure 5
gamma_df <- rbind(
  data.frame(AY = 1:m, gamma = as.numeric(dcl_det_par$inflat),
             model = "Baseline DCL"),
  data.frame(AY = 1:m, gamma = as.numeric(bdcl_det_par$inflat),
             model = "Baseline BDCL"),
  data.frame(AY = 1:m, gamma = as.numeric(dcl_test_par2$inflat.DCL),
             model = "DCL shocked"),
  data.frame(AY = 1:m, gamma = as.numeric(dcl_test_par2$inflat),
             model = "BDCL shocked")
)

gammaref <- ggplot(gamma_df, aes(AY, gamma, colour = model)) +
  geom_line() + geom_point(size = 2) +
  scale_colour_manual(values = c("Baseline DCL" = "grey60",
                                 "Baseline BDCL" = "grey30",
                                 "DCL shocked" = "skyblue3",
                                 "BDCL shocked" = "springgreen3")) +
  labs(title = expression("Severity Inflation " * gamma[i]), subtitle = sub,
       x = "Accident Year", y = expression(gamma[i])) + th
gammaref
################################################################################
#Appendix Parameter Plots for Baseline and Reference Scenario

Plot.dcl.par(dcl_test_par2, type.inflat = "BDCL")
Plot.clm.par(clm_sk)
Plot.dcl.par(testfit_bdcl_par, type.inflat = "BDCL")
Plot.clm.par(clm_base)


################################################################################


#Scenario Grid Plots


################################################################################
#Bias across scenarios, Figure 6
bias_df <- final %>%
  filter(metric == "RelTotalError") %>%
  filter(omega > 0) %>%
  mutate(
    model_group = ifelse(model == "BDCL", "BDCL", "CLM / DCL"),
    y_label = factor(paste0("y = ", y))
  )

#Use CLM results since CLM and DCL are almost identical
bias_plot <- bias_df %>%
  filter(!(model_group == "CLM / DCL" & model == "DCL"))

bias_sg <- ggplot(bias_plot, aes(omega, value, colour = y_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.75) +
  geom_point(size = 2) +
  facet_wrap(~model_group, nrow = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_colour_brewer(palette = "Set2", name = "Calendar-Year Shock Introduced") +
  labs(title = "Reserve Bias Across Shock Magnitudes",
       x = expression(omega ~ "(Shock Magnitude)"),
       y = "Reserve Bias") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )
bias_sg

################################################################################

#Worst CY Error scenario grid, Figure 7
worst_val <- final %>%
  filter(metric == "WorstCYError.WorstCYError") %>%
  select(y, omega, model, value)

worst_cy_num <- final %>%
  filter(metric == "WorstCYError.WorstCY") %>%
  select(y, omega, model, worst_cy = value)

worst_combined <- worst_val %>%
  left_join(worst_cy_num, by = c("y", "omega", "model")) %>%
  filter(omega > 0) %>%
  mutate(
    model = factor(model, levels = c("CLM", "DCL", "BDCL")),
    pct_value = value * 100,
    label = paste0(sprintf("%+.0f%%", pct_value), "\nCY ", round(worst_cy))
  )

max_abs <- max(abs(worst_combined$pct_value))
min_v <- min(worst_combined$pct_value)

ggplot(worst_combined, aes(x = factor(omega), y = factor(y), fill = pct_value)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 2.8, lineheight = 0.85) +
  facet_wrap(~model, ncol = 3) +
  scale_fill_gradient2(
    low = "steelblue",
    mid = "white",
    high = "brown2",
    midpoint = 0,
    limits = c(min_v, max_abs),
    name = "Worst CY\nError (%)"
  ) +
  labs(
    title = "Worst Calendar-Year Error by Scenario and Model",
    x = expression(omega ~ "(Shock Magnitude)"),
    y = "y  (Calendar-Year Shock is Introduced)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 12),

  )



################################################################################

#Wasserstein Distance Across Scenarios, Figure 8
wass_df_sg <- final %>%
  filter(metric == "WassersteinDistance") %>%
  filter(omega > 0) %>%
  mutate(
    model_group = ifelse(model == "BDCL", "BDCL", "CLM / DCL"),
    y_label = factor(paste0("y = ", y))
  ) %>%
  filter(!(model_group == "CLM / DCL" & model == "DCL"))

wass_sg <- ggplot(wass_df_sg, aes(omega, value, colour = y_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.75) +
  geom_point(size = 2) +
  facet_wrap(~model_group, nrow = 1) +
  scale_colour_brewer(palette = "Set2", name = "Calendar-Year Shock Introduced") +
  labs(title = "Cash Flow Timing Across Shock Magnitudes",
       x = expression(omega ~ "(Shock Magnitude)"),
       y = "1-Wasserstein Distance") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    strip.text       = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),

  )
wass_sg


################################################################################

#Gamma ratio plot across scenarios, Figure 9
gamma_list <- list()
for (i in 1:nrow(full_results)) {
  row <- full_results[i, ]
  if (row$omega == 0 || row$model == "CLM") next
  
  model_name <- as.character(row$model)
  base_row <- full_results %>%
    filter(model == model_name, omega == 0) %>%
    slice(1)
  
  if (nrow(base_row) == 0) next
  
  gamma_s <- as.numeric(if (model_name == "DCL") row$fitted[[1]]$inflat.DCL 
                        else row$fitted[[1]]$inflat)
  gamma_b <- as.numeric(if (model_name == "DCL") base_row$fitted[[1]]$inflat.DCL 
                        else base_row$fitted[[1]]$inflat)
  
  if (length(gamma_s) != m || length(gamma_b) != m) next
  
  gamma_list[[length(gamma_list) + 1]] <- data.frame(
    AY = 1:m, ratio = gamma_s / gamma_b,
    model = model_name, omega = as.numeric(row$omega), y = as.numeric(row$y)
  )
}

gamma_ratio_df <- bind_rows(gamma_list) %>%
  mutate(model = factor(model, levels = c("DCL", "BDCL")),
         omega_label = factor(paste0("ω = ", omega)))

ggplot(gamma_ratio_df, aes(AY, ratio - 1, colour = model, group = interaction(model, y))) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(alpha = 0.4, linewidth = 0.4) +
  stat_summary(aes(group = model), fun = median, geom = "line", linewidth = 1.1) +
  facet_wrap(~omega_label, nrow = 1) +
  scale_y_continuous(labels = scales::percent, breaks = scales::breaks_pretty(n = 10)) +
  scale_colour_manual(values = c(DCL = "skyblue3", BDCL = "springgreen3")) +
  labs(
    title = expression("Severity Parameter Response: " * gamma[i]^shocked / gamma[i]^baseline - 1),
    subtitle = "Bold = Median Across y Values; Light = Individual Scenarios",
    x = "Accident-Year",
    y = expression(gamma[i] ~ "Percentage Change (%)"),
    colour = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "grey80", fill = NA))

################################################################################

#CYRE Line Plot Across Scenarios, Figure 10
baseline_err <- cy_errors %>%
  filter(omega == 0) %>%
  mutate(model = factor(model, levels = c("CLM", "DCL", "BDCL")))

shocked_err <- cy_errors %>%
  filter(omega > 0) %>%
  mutate(
    scenario = paste0("y=", y, ", ω=", omega),
    model = factor(model, levels = c("CLM", "DCL", "BDCL"))
  )

p_cyerrorline <- ggplot() +
  #Shocked scenarios
  geom_line(
    data = shocked_err,
    aes(x = CY, y = value, group = scenario, alpha = omega),
    colour = "grey30", linewidth = 0.4
  ) +
  
  #Baseline
  geom_line(
    data = baseline_err,
    aes(x = CY, y = value),
    colour = "navy", linewidth = 1
  ) +
  
  #Reference line at 0
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_wrap(~model, nrow = 1) +
  scale_alpha_continuous(
    range = c(0.25, 1),
    breaks = c(0.5, 0.75, 1.0, 1.5),
    name = expression(omega)
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Calendar-Year Relative Error Across All Scenarios",
    subtitle = "Baseline (ω = 0) vs Shocked Scenarios",
    x = "Calendar-Year",
    y = "(pred - truth) / truth"
  ) +

  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

p_cyerrorline

################################################################################

wshortfall_results <- final %>% filter(metric == "WorstShortfall")

wshortfall_table <- wshortfall_results %>%
  select(y, omega, model, value) %>%
  pivot_wider(names_from = model, values_from = value) %>%
  mutate(across(c(CLM, DCL, BDCL), ~ sprintf("%.4f", .x)))
wshortfall_table
print(wshortfall_table, n = 28)
################################################################################

bias_results <- final %>% filter(metric == "RelTotalError")
bias_table <- bias_results %>%
  select(y, omega, model, value) %>%
  pivot_wider(names_from = model, values_from = value) %>%
  mutate(across(c(CLM, DCL, BDCL), ~ sprintf("%.4f", .x)))
bias_table
print(bias_table, n = 28)

