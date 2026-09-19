# ---- Header ------------------------------------------------------------
# Title: restoration_analysis_v05
# Author: HR
# Date: 2026-08-20
# Purpose: Analyse 2026 restoration-monitoring physicochemical data across
# six fixed sites and three sampling rounds.
#
# Input:
#   data/trt_clean_v06.rds
#
# Main outputs:
#   outputs/restoration_overall_summary.csv
#   outputs/restoration_round_summary.csv
#   outputs/restoration_site_summary.csv
#   outputs/restoration_sac_summary.csv
#   outputs/restoration_sac_by_round.csv
#   outputs/restoration_sac_by_site.csv
#   outputs/restoration_pairwise_rounds.csv
#   outputs/restoration_gis_site_summary.csv
#   outputs/restoration_gis_by_round.csv
#   outputs/appendix_B2_1_site_orthophosphate.csv
#   outputs/appendix_B2_2_pairwise_rounds.csv
#   outputs/appendix_B2_3_physchem_round_summary.csv
#   outputs/appendix_B2_4_physchem_observations.csv
#   outputs/figure_physchem_by_site.png
#
# Analytical notes:
# - Restoration monitoring is analysed independently from EA routine data.
# - Orthophosphate is expressed as mg/L as P.
# - Repeated measurements across the same six sites are compared using a
#   Friedman test, followed by paired Wilcoxon tests with BH adjustment.
# - SAC phosphorus targets are contextual annual-mean benchmarks, not formal
#   classifications of individual spot samples.

# ---- 0: Libraries and output directory --------------------------------

library(tidyverse)
library(lubridate)
library(patchwork)

dir.create(
  "outputs",
  showWarnings = FALSE,
  recursive = TRUE
)

theme_mease <- theme_minimal(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# ---- 1: Load and validate data -----------------------------------------

trt <- readRDS("data/trt_clean_v06.rds")

site_order <- c(
  "Upstream Gilwiskaw",
  "Wetland Gilwiskaw",
  "Upstream Mease",
  "Wetland Mease",
  "Confluence",
  "Downstream Mease"
)

site_code_order <- c("UG", "WG", "UM", "WM", "C", "DM")

trt_restoration <- trt %>%
  filter(
    source_programme == "WETLAND_RESTORATION"
  ) %>%
  mutate(
    site_name = factor(
      as.character(site_name),
      levels = site_order
    ),
    site_id_original = factor(
      site_id_original,
      levels = site_code_order
    ),
    survey_round = factor(
      survey_round,
      levels = sort(unique(survey_round)),
      labels = paste(
        "Round",
        sort(unique(survey_round))
      )
    )
  ) %>%
  arrange(
    date_sampled,
    site_id_original
  )

required_columns <- c(
  "site_id_original",
  "site_name",
  "tributary_name",
  "ngr",
  "easting",
  "northing",
  "survey_round",
  "date_sampled",
  "orthophosphate_as_p_mgl",
  "temperature_c",
  "dissolved_oxygen_mgl",
  "oxygen_saturation_percent",
  "ph"
)

stopifnot(
  all(required_columns %in% names(trt_restoration)),
  nrow(trt_restoration) == 18,
  n_distinct(trt_restoration$site_id_original) == 6,
  n_distinct(trt_restoration$survey_round) == 3,
  !any(is.na(trt_restoration$orthophosphate_as_p_mgl))
)

round_inventory <- trt_restoration %>%
  count(
    survey_round,
    date_sampled,
    name = "n_sites"
  )

duplicate_site_round_check <- trt_restoration %>%
  count(
    site_id_original,
    survey_round,
    name = "n"
  ) %>%
  filter(n != 1)

stopifnot(
  all(round_inventory$n_sites == 6),
  nrow(duplicate_site_round_check) == 0
)

# ---- 2: SAC phosphorus targets -----------------------------------------

trt_restoration <- trt_restoration %>%
  mutate(
    sac_p_target = case_when(
      tributary_name == "Gilwiskaw Brook" ~ 0.040,
      tributary_name == "River Mease" ~ 0.050,
      tributary_name == "Mease-Gilwiskaw Confluence" ~ 0.050,
      TRUE ~ NA_real_
    ),
    exceeds_sac_target =
      orthophosphate_as_p_mgl > sac_p_target,
    exceedance_mgl =
      orthophosphate_as_p_mgl - sac_p_target
  )

stopifnot(
  !any(is.na(trt_restoration$sac_p_target))
)

# ---- 3: Overall, round and site summaries ------------------------------

restoration_overall_summary <- trt_restoration %>%
  summarise(
    n_samples = n(),
    n_sites = n_distinct(site_id_original),
    n_rounds = n_distinct(survey_round),
    first_date = min(date_sampled),
    last_date = max(date_sampled),
    mean_orthophosphate =
      mean(orthophosphate_as_p_mgl),
    median_orthophosphate =
      median(orthophosphate_as_p_mgl),
    sd_orthophosphate =
      sd(orthophosphate_as_p_mgl),
    min_orthophosphate =
      min(orthophosphate_as_p_mgl),
    max_orthophosphate =
      max(orthophosphate_as_p_mgl)
  )

restoration_round_summary <- trt_restoration %>%
  group_by(
    survey_round,
    date_sampled
  ) %>%
  summarise(
    n_sites = n_distinct(site_id_original),
    mean_orthophosphate =
      mean(orthophosphate_as_p_mgl),
    median_orthophosphate =
      median(orthophosphate_as_p_mgl),
    sd_orthophosphate =
      sd(orthophosphate_as_p_mgl),
    min_orthophosphate =
      min(orthophosphate_as_p_mgl),
    max_orthophosphate =
      max(orthophosphate_as_p_mgl),
    mean_temperature =
      mean(temperature_c, na.rm = TRUE),
    mean_dissolved_oxygen =
      mean(dissolved_oxygen_mgl, na.rm = TRUE),
    mean_oxygen_saturation =
      mean(oxygen_saturation_percent, na.rm = TRUE),
    mean_ph =
      mean(ph, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(date_sampled)

restoration_site_summary <- trt_restoration %>%
  group_by(
    site_id_original,
    site_name,
    tributary_name,
    ngr,
    easting,
    northing
  ) %>%
  summarise(
    n_samples = n(),
    n_rounds = n_distinct(survey_round),
    mean_orthophosphate =
      mean(orthophosphate_as_p_mgl),
    median_orthophosphate =
      median(orthophosphate_as_p_mgl),
    sd_orthophosphate =
      sd(orthophosphate_as_p_mgl),
    min_orthophosphate =
      min(orthophosphate_as_p_mgl),
    max_orthophosphate =
      max(orthophosphate_as_p_mgl),
    mean_temperature =
      mean(temperature_c, na.rm = TRUE),
    mean_dissolved_oxygen =
      mean(dissolved_oxygen_mgl, na.rm = TRUE),
    mean_oxygen_saturation =
      mean(oxygen_saturation_percent, na.rm = TRUE),
    mean_ph =
      mean(ph, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(site_id_original)

# ---- 4: SAC target exceedance summaries --------------------------------

restoration_sac_summary <- trt_restoration %>%
  summarise(
    n_samples = n(),
    n_exceeding_target =
      sum(exceeds_sac_target),
    percent_exceeding_target =
      100 * mean(exceeds_sac_target)
  )

restoration_sac_by_round <- trt_restoration %>%
  group_by(
    survey_round,
    date_sampled
  ) %>%
  summarise(
    n_samples = n(),
    n_exceeding_target =
      sum(exceeds_sac_target),
    percent_exceeding_target =
      100 * mean(exceeds_sac_target),
    .groups = "drop"
  )

restoration_sac_by_site <- trt_restoration %>%
  group_by(
    site_id_original,
    site_name,
    tributary_name,
    sac_p_target
  ) %>%
  summarise(
    n_samples = n(),
    n_exceeding_target =
      sum(exceeds_sac_target),
    percent_exceeding_target =
      100 * mean(exceeds_sac_target),
    mean_exceedance_mgl =
      mean(exceedance_mgl),
    .groups = "drop"
  )

# ---- 5: Temporal inference across rounds -------------------------------

friedman_rounds <- friedman.test(
  orthophosphate_as_p_mgl ~
    survey_round |
    site_id_original,
  data = trt_restoration
)

pairwise_rounds <- pairwise.wilcox.test(
  x = trt_restoration$orthophosphate_as_p_mgl,
  g = trt_restoration$survey_round,
  paired = TRUE,
  p.adjust.method = "BH",
  exact = FALSE
)

pairwise_rounds_table <-
  as.data.frame(as.table(pairwise_rounds$p.value)) %>%
  as_tibble() %>%
  rename(
    round_1 = Var1,
    round_2 = Var2,
    bh_adjusted_p = Freq
  ) %>%
  filter(!is.na(bh_adjusted_p)) %>%
  arrange(bh_adjusted_p) %>%
  mutate(
    significant = bh_adjusted_p < 0.05
  )

# ---- 6: Physicochemical figure -----------------------------------------

physchem_plot_data <- trt_restoration %>%
  mutate(
    site_number = recode(
      as.character(site_id_original),
      UG = 1L,
      WG = 2L,
      UM = 3L,
      WM = 4L,
      C = 5L,
      DM = 6L
    ),
    site_number = factor(
      site_number,
      levels = 1:6
    )
  )

fig_temp <- ggplot(
  physchem_plot_data,
  aes(
    x = site_number,
    y = temperature_c,
    group = survey_round,
    shape = survey_round,
    linetype = survey_round
  )
) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2.5) +
  labs(
    x = NULL,
    y = "Water temperature (°C)",
    shape = "Monitoring round",
    linetype = "Monitoring round"
  ) +
  theme_mease

fig_do <- ggplot(
  physchem_plot_data,
  aes(
    x = site_number,
    y = dissolved_oxygen_mgl,
    group = survey_round,
    shape = survey_round,
    linetype = survey_round
  )
) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2.5) +
  labs(
    x = NULL,
    y = expression(
      paste("Dissolved oxygen (mg ", L^{-1}, ")")
    ),
    shape = "Monitoring round",
    linetype = "Monitoring round"
  ) +
  theme_mease

fig_sat <- ggplot(
  physchem_plot_data,
  aes(
    x = site_number,
    y = oxygen_saturation_percent,
    group = survey_round,
    shape = survey_round,
    linetype = survey_round
  )
) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2.5) +
  labs(
    x = "Site",
    y = "Oxygen saturation (%)",
    shape = "Monitoring round",
    linetype = "Monitoring round"
  ) +
  theme_mease

fig_ph <- ggplot(
  physchem_plot_data,
  aes(
    x = site_number,
    y = ph,
    group = survey_round,
    shape = survey_round,
    linetype = survey_round
  )
) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2.5) +
  labs(
    x = "Site",
    y = "pH",
    shape = "Monitoring round",
    linetype = "Monitoring round"
  ) +
  theme_mease

figure_physchem_by_site <- (
  fig_temp + fig_do +
    fig_sat + fig_ph
) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom"
  )

ggsave(
  "outputs/figure_physchem_by_site.png",
  figure_physchem_by_site,
  width = 9,
  height = 7,
  dpi = 300
)

# ---- 7: GIS-ready outputs ----------------------------------------------

restoration_gis_site <- restoration_site_summary %>%
  mutate(
    site_number = recode(
      as.character(site_id_original),
      UG = 1L,
      WG = 2L,
      UM = 3L,
      WM = 4L,
      C = 5L,
      DM = 6L
    ),
    sac_p_target = case_when(
      site_number %in% c(1, 2) ~ 0.040,
      site_number %in% c(3, 4, 5, 6) ~ 0.050
    ),
    target_status = if_else(
      median_orthophosphate > sac_p_target,
      "Target exceeded",
      "Target met"
    ),
    map_label = paste0(
      site_number,
      " (",
      sprintf("%.3f", median_orthophosphate),
      ")"
    )
  ) %>%
  arrange(site_number)

restoration_gis_round <- trt_restoration %>%
  transmute(
    site_id_original,
    site_name = as.character(site_name),
    tributary_name,
    ngr,
    easting,
    northing,
    survey_round = as.character(survey_round),
    date_sampled,
    orthophosphate_as_p_mgl,
    sac_p_target,
    exceeds_sac_target,
    exceedance_mgl,
    temperature_c,
    dissolved_oxygen_mgl,
    oxygen_saturation_percent,
    ph
  ) %>%
  arrange(
    date_sampled,
    site_id_original
  )

# ---- 8: Appendix B2 tables ---------------------------------------------

appendix_B2_1 <- restoration_site_summary %>%
  mutate(
    Site = recode(
      as.character(site_id_original),
      UG = 1L,
      WG = 2L,
      UM = 3L,
      WM = 4L,
      C = 5L,
      DM = 6L
    )
  ) %>%
  transmute(
    Site,
    `Site name` = as.character(site_name),
    n = n_samples,
    Mean = round(mean_orthophosphate, 3),
    Median = round(median_orthophosphate, 3),
    SD = round(sd_orthophosphate, 3),
    Minimum = round(min_orthophosphate, 3),
    Maximum = round(max_orthophosphate, 3)
  ) %>%
  arrange(Site)

appendix_B2_2 <- pairwise_rounds_table %>%
  transmute(
    Comparison = paste(
      as.character(round_1),
      as.character(round_2),
      sep = "–"
    ),
    `BH-adjusted p` = case_when(
      bh_adjusted_p < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", bh_adjusted_p)
    ),
    Significant = if_else(
      significant,
      "Yes",
      "No"
    )
  )

appendix_B2_3 <- restoration_round_summary %>%
  transmute(
    Round = as.character(survey_round),
    Date = date_sampled,
    `Mean temperature (°C)` =
      round(mean_temperature, 2),
    `Mean dissolved oxygen (mg L-1)` =
      round(mean_dissolved_oxygen, 2),
    `Mean oxygen saturation (%)` =
      round(mean_oxygen_saturation, 1),
    `Mean pH` =
      round(mean_ph, 2)
  )

appendix_B2_4 <- physchem_plot_data %>%
  transmute(
    Round = as.character(survey_round),
    Date = date_sampled,
    Site = as.integer(as.character(site_number)),
    `Temperature (°C)` = temperature_c,
    `Dissolved oxygen (mg L-1)` =
      dissolved_oxygen_mgl,
    `Oxygen saturation (%)` =
      oxygen_saturation_percent,
    pH = ph
  ) %>%
  arrange(Round, Site)

# ---- 9: Export ----------------------------------------------------------

write_csv(
  restoration_overall_summary,
  "outputs/restoration_overall_summary.csv"
)

write_csv(
  restoration_round_summary,
  "outputs/restoration_round_summary.csv"
)

write_csv(
  restoration_site_summary,
  "outputs/restoration_site_summary.csv"
)

write_csv(
  restoration_sac_summary,
  "outputs/restoration_sac_summary.csv"
)

write_csv(
  restoration_sac_by_round,
  "outputs/restoration_sac_by_round.csv"
)

write_csv(
  restoration_sac_by_site,
  "outputs/restoration_sac_by_site.csv"
)

write_csv(
  pairwise_rounds_table,
  "outputs/restoration_pairwise_rounds.csv"
)

write_csv(
  restoration_gis_site,
  "outputs/restoration_gis_site_summary.csv"
)

write_csv(
  restoration_gis_round,
  "outputs/restoration_gis_by_round.csv"
)

write_csv(
  appendix_B2_1,
  "outputs/appendix_B2_1_site_orthophosphate.csv"
)

write_csv(
  appendix_B2_2,
  "outputs/appendix_B2_2_pairwise_rounds.csv"
)

write_csv(
  appendix_B2_3,
  "outputs/appendix_B2_3_physchem_round_summary.csv"
)

write_csv(
  appendix_B2_4,
  "outputs/appendix_B2_4_physchem_observations.csv"
)

# ---- 10: Final QC / reporting summary ----------------------------------

cat("\nRestoration analysis complete\n")
cat("-----------------------------\n")

print(restoration_overall_summary)

cat("\nRound summary:\n")
print(restoration_round_summary)

cat("\nSite summary:\n")
print(restoration_site_summary)

cat("\nSAC target summary:\n")
print(restoration_sac_summary)

cat("\nFriedman test:\n")
print(friedman_rounds)

cat("\nPairwise paired Wilcoxon tests:\n")
print(pairwise_rounds_table)
