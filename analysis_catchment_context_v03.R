# ---- Header ------------------------------------------------------------
# Title: analysis_catchment_context_v03
# Author: HR
# Date: 2026-08-20
# Purpose: Analyse harmonised Environment Agency routine water-quality data
# to provide long-term catchment context for the River Mease restoration study.
# The workflow summarises spatial, temporal and seasonal variation in
# orthophosphate concentrations and produces the statistical, tabular, figure
# and GIS outputs used in the dissertation.
#
# Input:
#   data/trt_clean_v06.rds
#
# Main outputs:
#   outputs/catchment_site_summary.csv
#   outputs/catchment_tributary_summary.csv
#   outputs/catchment_year_summary.csv
#   outputs/catchment_season_summary.csv
#   outputs/catchment_pairwise_tributary.csv
#   outputs/catchment_pairwise_season.csv
#   outputs/site_gis_summary.csv
#   outputs/site_season_gis_summary.csv
#   outputs/appendix_B1_1.csv ... appendix_B1_5.csv
#   outputs/figure_3_2_tributary.png
#   outputs/figure_3_3_year.png
#   outputs/figure_3_4_season.png
#   outputs/figure_3_7_context_restoration.png
#
# Analytical scope:
#   - Historical EA routine monitoring: 2010-2025
#   - Contemporary restoration observations are retained only for Figure 3.7
#     and are not included in historical inferential analyses.
#
# Statistical interpretation:
#   Kruskal-Wallis and pairwise Wilcoxon tests are used as broad comparisons
#   of concentration distributions. Because some monitoring locations were
#   sampled repeatedly through time, these tests are not interpreted as fully
#   independent repeated-measures inference.

# ---- 0: Libraries and output directory --------------------------------

library(tidyverse)
library(lubridate)
library(rnrfa)

dir.create(
  "outputs",
  showWarnings = FALSE,
  recursive = TRUE
)

theme_mease <- theme_minimal(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold")
  )

# ---- 1: Load harmonised dataset ----------------------------------------

trt <- readRDS("data/trt_clean_v06.rds")

required_columns <- c(
  "site_id",
  "ngr",
  "tributary_name",
  "tributary_family",
  "date_sampled",
  "year",
  "source_programme",
  "orthophosphate_as_p_mgl"
)

stopifnot(
  all(required_columns %in% names(trt))
)

trt_context <- trt %>%
  filter(source_programme == "EA_ROUTINE")

trt_restoration <- trt %>%
  filter(source_programme == "WETLAND_RESTORATION")

# Confirm expected analytical structure
stopifnot(
  sum(!is.na(trt_context$orthophosphate_as_p_mgl)) == 633,
  nrow(trt_restoration) == 18,
  n_distinct(trt_restoration$survey_round) == 3
)

# ---- 2: Site-level spatial summaries -----------------------------------

site_summary <- trt_context %>%
  group_by(
    site_id,
    ngr,
    tributary_name,
    tributary_family
  ) %>%
  summarise(
    n = sum(!is.na(orthophosphate_as_p_mgl)),
    mean_ortho_p =
      mean(orthophosphate_as_p_mgl, na.rm = TRUE),
    median_ortho_p =
      median(orthophosphate_as_p_mgl, na.rm = TRUE),
    q25_ortho_p =
      quantile(orthophosphate_as_p_mgl, 0.25, na.rm = TRUE),
    q75_ortho_p =
      quantile(orthophosphate_as_p_mgl, 0.75, na.rm = TRUE),
    max_ortho_p =
      max(orthophosphate_as_p_mgl, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  site_summary,
  "outputs/catchment_site_summary.csv"
)

# ---- 3: QGIS-ready site summaries --------------------------------------

# Convert National Grid References to British National Grid Easting/Northing.
# These exports contain one row per monitoring site.

site_gis_export <- site_summary %>%
  rowwise() %>%
  mutate(
    coords = list(osg_parse(ngr)),
    easting = coords$easting,
    northing = coords$northing
  ) %>%
  ungroup() %>%
  select(-coords) %>%
  arrange(desc(median_ortho_p))

write_csv(
  site_gis_export,
  "outputs/site_gis_summary.csv"
)

# Seasonal site-level export for QGIS
site_season_gis_export <- trt_context %>%
  mutate(
    month = month(date_sampled),
    season = case_when(
      month %in% c(12, 1, 2) ~ "winter",
      month %in% c(3, 4, 5) ~ "spring",
      month %in% c(6, 7, 8) ~ "summer",
      month %in% c(9, 10, 11) ~ "autumn",
      TRUE ~ NA_character_
    ),
    season = factor(
      season,
      levels = c("spring", "summer", "autumn", "winter")
    )
  ) %>%
  filter(
    !is.na(season),
    !is.na(orthophosphate_as_p_mgl)
  ) %>%
  group_by(
    site_id,
    ngr,
    tributary_name,
    tributary_family,
    season
  ) %>%
  summarise(
    n_obs = n(),
    mean_ortho =
      mean(orthophosphate_as_p_mgl),
    median_ortho =
      median(orthophosphate_as_p_mgl),
    max_ortho =
      max(orthophosphate_as_p_mgl),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    coords = list(osg_parse(ngr)),
    easting = coords$easting,
    northing = coords$northing
  ) %>%
  ungroup() %>%
  select(-coords) %>%
  arrange(site_id, season)

write_csv(
  site_season_gis_export,
  "outputs/site_season_gis_summary.csv"
)

# ---- 4: Tributary-family summaries -------------------------------------

tributary_summary <- trt_context %>%
  filter(
    !is.na(tributary_family),
    !is.na(orthophosphate_as_p_mgl)
  ) %>%
  group_by(tributary_family) %>%
  summarise(
    n = n(),
    n_years = n_distinct(year),
    mean_ortho_p =
      mean(orthophosphate_as_p_mgl),
    median_ortho_p =
      median(orthophosphate_as_p_mgl),
    max_ortho_p =
      max(orthophosphate_as_p_mgl),
    .groups = "drop"
  ) %>%
  arrange(median_ortho_p)

write_csv(
  tributary_summary,
  "outputs/catchment_tributary_summary.csv"
)

# Restrict inferential comparisons to tributary families represented in at
# least three monitoring years, reducing the influence of very sparse groups.
tributaries_multi_year <- tributary_summary %>%
  filter(n_years >= 3) %>%
  pull(tributary_family)

trt_multi_year <- trt_context %>%
  filter(
    tributary_family %in% tributaries_multi_year,
    !is.na(orthophosphate_as_p_mgl)
  )

stopifnot(
  n_distinct(trt_multi_year$tributary_family) >= 3
)

# ---- 5: Tributary-family tests -----------------------------------------

tributary_kruskal <- kruskal.test(
  orthophosphate_as_p_mgl ~ tributary_family,
  data = trt_multi_year
)

tributary_pairwise <- pairwise.wilcox.test(
  x = trt_multi_year$orthophosphate_as_p_mgl,
  g = trt_multi_year$tributary_family,
  p.adjust.method = "BH",
  exact = FALSE
)

tributary_pairwise_table <-
  as.data.frame(as.table(tributary_pairwise$p.value)) %>%
  as_tibble() %>%
  rename(
    tributary_1 = Var1,
    tributary_2 = Var2,
    bh_adjusted_p = Freq
  ) %>%
  filter(!is.na(bh_adjusted_p)) %>%
  arrange(bh_adjusted_p) %>%
  mutate(
    significant = bh_adjusted_p < 0.05
  )

write_csv(
  tributary_pairwise_table,
  "outputs/catchment_pairwise_tributary.csv"
)

# ---- 6: Figure 3.2 - tributary variation -------------------------------

tributary_labels <- c(
  "MEASE - UPPER MEASE"       = "Mease - Upper Mease (1)",
  "ASHBY CANAL"               = "Ashby Canal (2)",
  "MEASE - SNARESTONE STW"    = "Mease - Snarestone STW (3)",
  "GILWISKAW BROOK"           = "Gilwiskaw Brook (4)",
  "APPLEBY BROOK"             = "Appleby Brook (5)",
  "MEASE - BIRDSHILL ROAD"    = "Mease - Birdshill Road (6)",
  "SHELL BROOK"               = "Shell Brook (7)",
  "SALTERSFORD BROOK"         = "Saltersford Brook (8)",
  "MEASE - STRETTON BRIDGE"   = "Mease - Stretton Bridge (9)",
  "HOOBOROUGH BROOK"          = "Hooborough Brook (10)",
  "KES BROOK"                 = "Kes Brook (11)",
  "CHILCOTE BROOK"            = "Chilcote Brook (12)",
  "MEASE - CLIFTON CAMPVILLE" = "Mease - Clifton Campville (13)",
  "SEAL BROOK"                = "Seal Brook (14)",
  "WEST BROOK"                = "West Brook (15)",
  "HARLASTON BROOK"           = "Harlaston Brook (16)",
  "PESSAL BROOK"              = "Pessal Brook (17)",
  "MEASE - CROXALL"           = "Mease - Croxall (18)",
  "MEASE - MAIN STEM"         = "Mease - Main Stem (19)"
)

fig_3_2 <- trt_context %>%
  filter(
    !is.na(tributary_family),
    !is.na(orthophosphate_as_p_mgl)
  ) %>%
  ggplot(
    aes(
      x = reorder(
        tributary_family,
        orthophosphate_as_p_mgl,
        median,
        na.rm = TRUE
      ),
      y = orthophosphate_as_p_mgl
    )
  ) +
  geom_boxplot(outlier.shape = NA) +
  coord_flip() +
  scale_x_discrete(labels = tributary_labels) +
  scale_y_log10() +
  labs(
    x = "Tributary family",
    y = expression(
      paste(
        "Orthophosphate (mg ", L^{-1},
        " as P; log"[10], " scale)"
      )
    )
  ) +
  theme_mease

ggsave(
  "outputs/figure_3_2_tributary.png",
  fig_3_2,
  width = 8,
  height = 7,
  dpi = 300
)

# ---- 7: Annual temporal summaries --------------------------------------

year_summary <- trt_context %>%
  filter(!is.na(orthophosphate_as_p_mgl)) %>%
  group_by(year) %>%
  summarise(
    n = n(),
    mean_ortho_p =
      mean(orthophosphate_as_p_mgl),
    median_ortho_p =
      median(orthophosphate_as_p_mgl),
    .groups = "drop"
  )

write_csv(
  year_summary,
  "outputs/catchment_year_summary.csv"
)

# Simple raw-observation model.
# Interpreted as an exploratory directional-trend model because monitoring
# effort is uneven and observations include repeated sampling at sites.
trend_model_raw <- lm(
  orthophosphate_as_p_mgl ~ year,
  data = trt_context
)

# Supplementary model using tributary-family annual means to reduce the
# dominance of heavily sampled site/year combinations.
trt_year_family <- trt_multi_year %>%
  group_by(
    tributary_family,
    year
  ) %>%
  summarise(
    mean_ortho =
      mean(orthophosphate_as_p_mgl),
    .groups = "drop"
  )

trend_model_family <- lm(
  log(mean_ortho + 0.001) ~ year,
  data = trt_year_family
)

# Compact model summaries for transparent reporting
trend_summary <- tibble(
  model = c(
    "individual_observations",
    "tributary_family_annual_means"
  ),
  n = c(
    nobs(trend_model_raw),
    nobs(trend_model_family)
  ),
  slope = c(
    unname(coef(trend_model_raw)["year"]),
    unname(coef(trend_model_family)["year"])
  ),
  p_value = c(
    summary(trend_model_raw)$coefficients["year", "Pr(>|t|)"],
    summary(trend_model_family)$coefficients["year", "Pr(>|t|)"]
  ),
  r_squared = c(
    summary(trend_model_raw)$r.squared,
    summary(trend_model_family)$r.squared
  )
)

write_csv(
  trend_summary,
  "outputs/catchment_temporal_trend_models.csv"
)

# ---- 8: Figure 3.3 - interannual variation -----------------------------

year_labels <- year_summary %>%
  mutate(
    label = paste0(year, "\n(n = ", n, ")")
  ) %>%
  select(year, label) %>%
  deframe()

fig_3_3 <- trt_context %>%
  filter(!is.na(orthophosphate_as_p_mgl)) %>%
  ggplot(
    aes(
      x = factor(year),
      y = orthophosphate_as_p_mgl
    )
  ) +
  geom_boxplot(outlier.shape = NA) +
  scale_x_discrete(labels = year_labels) +
  scale_y_log10() +
  labs(
    x = "Year",
    y = expression(
      paste(
        "Orthophosphate (mg ", L^{-1},
        " as P; log"[10], " scale)"
      )
    )
  ) +
  theme_mease

ggsave(
  "outputs/figure_3_3_year.png",
  fig_3_3,
  width = 8,
  height = 5,
  dpi = 300
)

# ---- 9: Seasonal summaries and tests -----------------------------------

trt_context_seasonal <- trt_context %>%
  mutate(
    month = month(date_sampled),
    season = case_when(
      month %in% c(12, 1, 2) ~ "winter",
      month %in% c(3, 4, 5) ~ "spring",
      month %in% c(6, 7, 8) ~ "summer",
      month %in% c(9, 10, 11) ~ "autumn",
      TRUE ~ NA_character_
    ),
    season = factor(
      season,
      levels = c(
        "spring",
        "summer",
        "autumn",
        "winter"
      )
    )
  ) %>%
  filter(
    !is.na(season),
    !is.na(orthophosphate_as_p_mgl)
  )

season_summary <- trt_context_seasonal %>%
  group_by(season) %>%
  summarise(
    n = n(),
    mean_ortho =
      mean(orthophosphate_as_p_mgl),
    median_ortho =
      median(orthophosphate_as_p_mgl),
    sd_ortho =
      sd(orthophosphate_as_p_mgl),
    min_ortho =
      min(orthophosphate_as_p_mgl),
    max_ortho =
      max(orthophosphate_as_p_mgl),
    .groups = "drop"
  )

write_csv(
  season_summary,
  "outputs/catchment_season_summary.csv"
)

season_kruskal <- kruskal.test(
  orthophosphate_as_p_mgl ~ season,
  data = trt_context_seasonal
)

season_pairwise <- pairwise.wilcox.test(
  x = trt_context_seasonal$orthophosphate_as_p_mgl,
  g = trt_context_seasonal$season,
  p.adjust.method = "BH",
  exact = FALSE
)

season_pairwise_table <-
  as.data.frame(as.table(season_pairwise$p.value)) %>%
  as_tibble() %>%
  rename(
    season_1 = Var1,
    season_2 = Var2,
    bh_adjusted_p = Freq
  ) %>%
  filter(!is.na(bh_adjusted_p)) %>%
  arrange(bh_adjusted_p) %>%
  mutate(
    significant = bh_adjusted_p < 0.05
  )

write_csv(
  season_pairwise_table,
  "outputs/catchment_pairwise_season.csv"
)

# ---- 10: Figure 3.4 - seasonal variation -------------------------------

season_labels <- season_summary %>%
  mutate(
    label = paste0(
      str_to_title(as.character(season)),
      "\n(n = ", n, ")"
    )
  ) %>%
  select(season, label) %>%
  deframe()

fig_3_4 <- ggplot(
  trt_context_seasonal,
  aes(
    x = season,
    y = orthophosphate_as_p_mgl
  )
) +
  geom_boxplot() +
  scale_x_discrete(labels = season_labels) +
  scale_y_log10() +
  labs(
    x = "Season",
    y = expression(
      paste(
        "Orthophosphate (mg ", L^{-1},
        " as P; log"[10], " scale)"
      )
    )
  ) +
  theme_mease

ggsave(
  "outputs/figure_3_4_season.png",
  fig_3_4,
  width = 7,
  height = 5,
  dpi = 300
)

# ---- 11: Figure 3.7 - historical context vs 2026 -----------------------

historic_annual <- trt_context %>%
  filter(!is.na(orthophosphate_as_p_mgl)) %>%
  transmute(
    plot_group = as.character(year),
    orthophosphate_as_p_mgl
  )

historic_summer <- trt_context %>%
  filter(
    month(date_sampled) %in% c(6, 7, 8),
    !is.na(orthophosphate_as_p_mgl)
  ) %>%
  transmute(
    plot_group = "2010-2025 Summer",
    orthophosphate_as_p_mgl
  )

restoration_2026 <- trt_restoration %>%
  filter(!is.na(orthophosphate_as_p_mgl)) %>%
  transmute(
    plot_group = "2026 Restoration Summer",
    orthophosphate_as_p_mgl
  )

catchment_restoration_comparison <- bind_rows(
  historic_annual,
  historic_summer,
  restoration_2026
) %>%
  mutate(
    plot_group = factor(
      plot_group,
      levels = c(
        "2010",
        "2020",
        "2021",
        "2022",
        "2025",
        "2010-2025 Summer",
        "2026 Restoration Summer"
      )
    )
  )

comparison_n <- catchment_restoration_comparison %>%
  count(plot_group, .drop = FALSE)

comparison_labels <- setNames(
  paste0(
    comparison_n$plot_group,
    "\n(n = ", comparison_n$n, ")"
  ),
  comparison_n$plot_group
)

comparison_labels["2010-2025 Summer"] <- paste0(
  "2010-2025\nSummer\n(n = ",
  comparison_n$n[
    comparison_n$plot_group == "2010-2025 Summer"
  ],
  ")"
)

comparison_labels["2026 Restoration Summer"] <- paste0(
  "2026 Restoration\nSummer\n(n = ",
  comparison_n$n[
    comparison_n$plot_group == "2026 Restoration Summer"
  ],
  ")"
)

fig_3_7 <- ggplot(
  catchment_restoration_comparison,
  aes(
    x = plot_group,
    y = orthophosphate_as_p_mgl
  )
) +
  geom_boxplot(outlier.shape = NA) +
  geom_vline(
    xintercept = 6.5,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  scale_x_discrete(labels = comparison_labels) +
  scale_y_log10() +
  labs(
    x = NULL,
    y = expression(
      paste(
        "Orthophosphate (mg ", L^{-1},
        " as P; log"[10], " scale)"
      )
    )
  ) +
  theme_mease

ggsave(
  "outputs/figure_3_7_context_restoration.png",
  fig_3_7,
  width = 9,
  height = 5,
  dpi = 300
)

# ---- 12: Historical WQ periods at focal locations ----------------------

# Retain the historical site-period summary used to document temporal coverage
# at focal locations corresponding to the contemporary monitoring network.

historic_wq_periods <- trt_context %>%
  filter(
    site_id %in% c(
      4,           # Birds Hill
      43, 19,      # B-road Gilwiskaw
      45,          # Upper Mease
      31, 30, 29,  # Confluence
      5            # Restoration Reach Mease
    )
  ) %>%
  group_by(
    site_id,
    ngr,
    tributary_name,
    tributary_family
  ) %>%
  summarise(
    first_year = min(year, na.rm = TRUE),
    last_year = max(year, na.rm = TRUE),
    n_years_sampled = n_distinct(year),
    n_obs = sum(!is.na(orthophosphate_as_p_mgl)),
    mean_ortho =
      mean(orthophosphate_as_p_mgl, na.rm = TRUE),
    median_ortho =
      median(orthophosphate_as_p_mgl, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    historic_wq_period =
      paste0(first_year, "-", last_year)
  ) %>%
  arrange(site_id, tributary_name)

write_csv(
  historic_wq_periods,
  "outputs/historic_wq_periods_focal_sites.csv"
)

# ---- 13: Appendix B1 tables --------------------------------------------

# B1.1 - tributary-family descriptive statistics
appendix_B1_1 <- tributary_summary %>%
  transmute(
    `Tributary family` = tributary_family,
    `n` = n,
    `Monitoring years` = n_years,
    `Mean` = round(mean_ortho_p, 3),
    `Median` = round(median_ortho_p, 3),
    `Maximum` = round(max_ortho_p, 3)
  )

write_csv(
  appendix_B1_1,
  "outputs/appendix_B1_1.csv"
)

# B1.2 - pairwise tributary-family comparisons
appendix_B1_2 <- tributary_pairwise_table %>%
  transmute(
    `Tributary family 1` = tributary_1,
    `Tributary family 2` = tributary_2,
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

write_csv(
  appendix_B1_2,
  "outputs/appendix_B1_2.csv"
)

# B1.3 - annual descriptive statistics
appendix_B1_3 <- year_summary %>%
  transmute(
    Year = year,
    `n` = n,
    Mean = round(mean_ortho_p, 3),
    Median = round(median_ortho_p, 3)
  )

write_csv(
  appendix_B1_3,
  "outputs/appendix_B1_3.csv"
)

# B1.4 - seasonal descriptive statistics
appendix_B1_4 <- season_summary %>%
  transmute(
    Season = str_to_title(as.character(season)),
    `n` = n,
    Mean = round(mean_ortho, 3),
    Median = round(median_ortho, 3),
    SD = round(sd_ortho, 3),
    Minimum = round(min_ortho, 3),
    Maximum = round(max_ortho, 3)
  )

write_csv(
  appendix_B1_4,
  "outputs/appendix_B1_4.csv"
)

# B1.5 - pairwise seasonal comparisons
appendix_B1_5 <- season_pairwise_table %>%
  transmute(
    `Season 1` = str_to_title(as.character(season_1)),
    `Season 2` = str_to_title(as.character(season_2)),
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

write_csv(
  appendix_B1_5,
  "outputs/appendix_B1_5.csv"
)

# ---- 14: Final analytical QC summary -----------------------------------

cat("\nCatchment-context analysis complete\n")
cat("-----------------------------------\n")
cat(
  "Historic usable orthophosphate observations:",
  sum(!is.na(trt_context$orthophosphate_as_p_mgl)),
  "\n"
)
cat(
  "Historic monitoring sites:",
  n_distinct(trt_context$site_id),
  "\n"
)
cat(
  "Tributary families:",
  n_distinct(
    trt_context$tributary_family[
      !is.na(trt_context$tributary_family)
    ]
  ),
  "\n"
)
cat(
  "Multi-year tributary families used inferentially:",
  length(tributaries_multi_year),
  "\n\n"
)

cat("Tributary Kruskal-Wallis test:\n")
print(tributary_kruskal)

cat("\nSeasonal Kruskal-Wallis test:\n")
print(season_kruskal)

cat("\nTemporal trend models:\n")
print(trend_summary)

cat("\nAnnual summary:\n")
print(year_summary)

cat("\nSeasonal summary:\n")
print(season_summary)

cat("\nSignificant tributary pairwise comparisons:\n")
print(
  tributary_pairwise_table %>%
    filter(significant)
)

cat("\nSeasonal pairwise comparisons:\n")
print(season_pairwise_table)
