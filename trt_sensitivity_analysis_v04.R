# ---- Header ------------------------------------------------------------
# Title: trt_sensitivity_analysis_v04
# Author: HR
# Date: 2026-08-20
# Purpose: Test whether the principal catchment-scale phosphorus conclusions
# are sensitive to alternative unit assumptions for orthophosphate in the
# 2021 and 2025 monitoring datasets.
#
# Input:
#   data/trt_clean_v06.rds
#
# Main outputs:
#   outputs/sensitivity_year_compare.csv
#   outputs/sensitivity_tributary_ranks.csv
#   outputs/sensitivity_site_rank_correlations.csv
#   outputs/sensitivity_trend_summary.csv
#   outputs/sensitivity_top5_overlap.csv
#   outputs/appendix_B1_6.csv
#   outputs/appendix_B1_7.csv
#
# Scenario definitions:
#   A = Current harmonisation workflow
#   B = 2021 values treated as already reported as P
#   C = 2025 values treated as already reported as P
#   D = 2021 and 2025 both treated as already reported as P
#
# Interpretation:
#   The analysis is diagnostic. It does not overwrite the cleaned dataset.
#   The 2021 reporting basis remains uncertain; 2025 is explicitly reported as
#   PO4 in the source metadata, so scenarios C and D are deliberately
#   counterfactual robustness checks.

# ---- 0: Libraries and output directory --------------------------------

library(tidyverse)

dir.create(
  "outputs",
  showWarnings = FALSE,
  recursive = TRUE
)

# ---- 1: Load historical EA dataset -------------------------------------

trt_base <- readRDS(
  "data/trt_clean_v06.rds"
) %>%
  filter(
    source_programme == "EA_ROUTINE"
  )

stopifnot(
  sum(!is.na(trt_base$orthophosphate_as_p_mgl)) == 633
)

# The original PO4-scale values are retained in
# orthophosphate_as_po4_mgl for 2021 and 2025, allowing alternative
# reporting assumptions to be reconstructed without editing the source data.

# ---- 2: Create sensitivity scenarios -----------------------------------

trt_A <- trt_base %>%
  mutate(
    scenario = "A_current_pipeline"
  )

trt_B <- trt_base %>%
  mutate(
    orthophosphate_as_p_mgl = case_when(
      year == 2021 &
        !is.na(orthophosphate_as_po4_mgl) ~
        orthophosphate_as_po4_mgl,
      TRUE ~ orthophosphate_as_p_mgl
    ),
    scenario = "B_2021_no_conversion"
  )

trt_C <- trt_base %>%
  mutate(
    orthophosphate_as_p_mgl = case_when(
      year == 2025 &
        !is.na(orthophosphate_as_po4_mgl) ~
        orthophosphate_as_po4_mgl,
      TRUE ~ orthophosphate_as_p_mgl
    ),
    scenario = "C_2025_no_conversion"
  )

trt_D <- trt_base %>%
  mutate(
    orthophosphate_as_p_mgl = case_when(
      year %in% c(2021, 2025) &
        !is.na(orthophosphate_as_po4_mgl) ~
        orthophosphate_as_po4_mgl,
      TRUE ~ orthophosphate_as_p_mgl
    ),
    scenario = "D_2021_2025_no_conversion"
  )

trt_sens <- bind_rows(
  trt_A,
  trt_B,
  trt_C,
  trt_D
)

# Confirm four complete scenario copies were created
stopifnot(
  nrow(trt_sens) == 4 * nrow(trt_base),
  n_distinct(trt_sens$scenario) == 4
)

# ---- 3: Annual summaries -----------------------------------------------

year_compare <- trt_sens %>%
  filter(
    !is.na(orthophosphate_as_p_mgl)
  ) %>%
  group_by(
    scenario,
    year
  ) %>%
  summarise(
    n = n(),
    mean_ortho =
      mean(orthophosphate_as_p_mgl),
    median_ortho =
      median(orthophosphate_as_p_mgl),
    .groups = "drop"
  )

write_csv(
  year_compare,
  "outputs/sensitivity_year_compare.csv"
)

# ---- 4: Tributary-family rankings --------------------------------------

trib_compare <- trt_sens %>%
  filter(
    !is.na(tributary_family),
    !is.na(orthophosphate_as_p_mgl)
  ) %>%
  group_by(
    scenario,
    tributary_family
  ) %>%
  summarise(
    median_ortho =
      median(orthophosphate_as_p_mgl),
    .groups = "drop"
  ) %>%
  group_by(scenario) %>%
  arrange(
    desc(median_ortho),
    .by_group = TRUE
  ) %>%
  mutate(
    rank = row_number()
  ) %>%
  ungroup()

write_csv(
  trib_compare,
  "outputs/sensitivity_tributary_ranks.csv"
)

# ---- 5: Site-level rank stability --------------------------------------

site_ranks <- trt_sens %>%
  filter(
    !is.na(orthophosphate_as_p_mgl)
  ) %>%
  group_by(
    scenario,
    site_id
  ) %>%
  summarise(
    median_ortho =
      median(orthophosphate_as_p_mgl),
    .groups = "drop"
  )

site_wide <- site_ranks %>%
  pivot_wider(
    names_from = scenario,
    values_from = median_ortho
  )

cor_matrix <- cor(
  site_wide %>%
    select(-site_id),
  use = "pairwise.complete.obs",
  method = "spearman"
)

site_rank_correlations <- tibble(
  scenario = c(
    "A_current_pipeline",
    "B_2021_no_conversion",
    "C_2025_no_conversion",
    "D_2021_2025_no_conversion"
  ),
  spearman_rho_vs_A = c(
    1,
    unname(
      cor_matrix[
        "A_current_pipeline",
        "B_2021_no_conversion"
      ]
    ),
    unname(
      cor_matrix[
        "A_current_pipeline",
        "C_2025_no_conversion"
      ]
    ),
    unname(
      cor_matrix[
        "A_current_pipeline",
        "D_2021_2025_no_conversion"
      ]
    )
  )
)

write_csv(
  site_rank_correlations,
  "outputs/sensitivity_site_rank_correlations.csv"
)

# ---- 6: Temporal trend stability ---------------------------------------

trend_summary <- trt_sens %>%
  filter(
    !is.na(orthophosphate_as_p_mgl)
  ) %>%
  group_by(scenario) %>%
  group_modify(
    ~ {
      model <- lm(
        orthophosphate_as_p_mgl ~ year,
        data = .x
      )

      tibble(
        n = nobs(model),
        slope =
          unname(coef(model)["year"]),
        p_value =
          summary(model)$coefficients[
            "year",
            "Pr(>|t|)"
          ],
        r_squared =
          summary(model)$r.squared
      )
    }
  ) %>%
  ungroup()

write_csv(
  trend_summary,
  "outputs/sensitivity_trend_summary.csv"
)

# ---- 7: Top-five tributary overlap -------------------------------------

top5_by_scenario <- trib_compare %>%
  group_by(scenario) %>%
  slice_min(
    order_by = rank,
    n = 5,
    with_ties = FALSE
  ) %>%
  summarise(
    top5 = list(tributary_family),
    .groups = "drop"
  )

reference_top5 <- top5_by_scenario %>%
  filter(
    scenario == "A_current_pipeline"
  ) %>%
  pull(top5) %>%
  .[[1]]

top5_overlap <- top5_by_scenario %>%
  rowwise() %>%
  mutate(
    n_overlap_with_A =
      length(intersect(top5, reference_top5)),
    overlap_fraction =
      paste0(n_overlap_with_A, "/5"),
    top5_text =
      paste(top5, collapse = " | ")
  ) %>%
  ungroup() %>%
  select(
    scenario,
    n_overlap_with_A,
    overlap_fraction,
    top5_text
  )

write_csv(
  top5_overlap,
  "outputs/sensitivity_top5_overlap.csv"
)

# ---- 8: Appendix B1.6 --------------------------------------------------

scenario_labels <- c(
  "A_current_pipeline" =
    "A: Current pipeline",
  "B_2021_no_conversion" =
    "B: 2021 not converted",
  "C_2025_no_conversion" =
    "C: 2025 not converted",
  "D_2021_2025_no_conversion" =
    "D: 2021 and 2025 not converted"
)

appendix_B1_6 <- year_compare %>%
  mutate(
    Scenario = recode(
      scenario,
      !!!scenario_labels
    )
  ) %>%
  transmute(
    Scenario,
    Year = year,
    `n` = n,
    Mean = round(mean_ortho, 3),
    Median = round(median_ortho, 3)
  )

write_csv(
  appendix_B1_6,
  "outputs/appendix_B1_6.csv"
)

# ---- 9: Appendix B1.7 --------------------------------------------------

# Generate robustness statistics directly from the analytical objects above.
appendix_B1_7 <- trend_summary %>%
  left_join(
    site_rank_correlations,
    by = "scenario"
  ) %>%
  left_join(
    top5_overlap %>%
      select(
        scenario,
        overlap_fraction
      ),
    by = "scenario"
  ) %>%
  mutate(
    Scenario = recode(
      scenario,
      !!!scenario_labels
    )
  ) %>%
  transmute(
    Scenario,
    `Site-rank rho` =
      round(spearman_rho_vs_A, 3),
    `Top-five overlap` =
      overlap_fraction,
    Slope =
      round(slope, 5),
    `p` =
      case_when(
        p_value < 0.001 ~ "<0.001",
        TRUE ~ sprintf("%.3f", p_value)
      ),
    `R-squared` =
      round(r_squared, 3)
  )

write_csv(
  appendix_B1_7,
  "outputs/appendix_B1_7.csv"
)

# ---- 10: Final QC / reporting summary ----------------------------------

cat("\nSensitivity analysis complete\n")
cat("-----------------------------\n")

cat("\nSite-rank correlations with Scenario A:\n")
print(site_rank_correlations)

cat("\nTemporal trend sensitivity:\n")
print(trend_summary)

cat("\nTop-five tributary overlap with Scenario A:\n")
print(
  top5_overlap %>%
    select(
      scenario,
      overlap_fraction
    )
)

cat("\nAppendix B1.7:\n")
print(appendix_B1_7)
