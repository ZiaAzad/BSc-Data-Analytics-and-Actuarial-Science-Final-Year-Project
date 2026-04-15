library(DCL)
library(tidyverse)
library(magrittr)
library(emdist)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)

source("~/Documents/DISSERTATION/Project/Data and Utils.R")
source("~/Documents/DISSERTATION/Project/Deterministic.R")
source("~/Documents/DISSERTATION/Project/Shock.R")
source("~/Documents/DISSERTATION/Project/Metrics.R")
source("~/Documents/DISSERTATION/Project/Scenario.R")


results <- vector("list", nrow(scenario_matrix))
for (i in seq_len(nrow(scenario_matrix))){
  results[[i]] <- run_scenario(scenario_matrix$y[i], scenario_matrix$omega[i])
}
full_results <- bind_rows(results)

final <- full_results %>%
  mutate(
    metrics_df = map(metrics, ~ tibble(
      metric = names(.x),
      value  = as.numeric(.x)
    ))
  ) %>%
  select(y, omega, model, metrics_df) %>%
  unnest(metrics_df) %>%
  mutate(
    value = ifelse(!is.finite(value), NA_real_, ifelse(abs(value) < tol, 0, as.numeric(value)))
  )

print(final, n = 441)


cy_errors <- full_results %>%
  transmute(y, omega, model, cy_rel) %>%
  mutate(cy_rel = map(cy_rel, ~ as.data.frame(.x))) %>%
  unnest(cy_rel) %>%
  rename(CY = CY, value = value) %>%
  mutate(CY = as.numeric(CY)) %>%
  mutate(
    value = ifelse(!is.finite(value), NA_real_, ifelse(abs(value) < tol, 0, value))
  )

cy_errors
print(cy_errors, n = 1134)
