# ============================================================================
# River Mease floodplain restoration | Historical macroinvertebrates
# Script 03: Sampling effort and taxon-richness comparisons
# File: macro_analysis_step02_v02.R
# ============================================================================
# PURPOSE
#   Summarise monitoring effort and per-sample taxon richness across the full
#   historical catchment dataset, six study sites and three long-term sites.
#   Produce descriptive plots and exploratory, unadjusted-for-year site tests.
#
# RUN ORDER
#   1. macro_historical_v03.R
#   2. macro_analysis_step01_v02.R
#   3. This script, from the EA_data project root (not its scripts folder).
#
# INPUTS
#   data/processed/macro_metrics_v01.rds
#   data/processed/macro_study_metrics_v01.rds
#
# OUTPUTS
#   Original step02 CSV/PNG filenames are retained for downstream compatibility.
#   data/processed/macro_longterm_v01.rds
#   Additional step02 test summaries and site table exported as CSV.
#
# ANALYTICAL DECISIONS / LIMITATIONS
#   Richness is the count of source taxon columns with recorded values > 0,
#   calculated in Script 02. Mixed historic abundance scales prohibit treating
#   summed abundance as a comparable ecological response across years.
#   Sampling effort varies by site and year; site tests are exploratory
#   unadjusted comparisons, not estimates of restoration effects.
#   The 2021 period boundary is descriptive, not a causal before/after design.
#   Historical and contemporary sites are not necessarily geographically
#   identical. A line between annual medians is visual, not a fitted trend.
#   Site labels follow the six-site lookup established in Script 02.
#
# STATUS
#   Revised from macro_analysis_step02_v01.R. Run and compare outputs against
#   the submitted analysis before treating this version as validated.
# ============================================================================

# ---- 0. Packages, paths and input checks ---------------------------------
library(tidyverse)
library(knitr)

processed_dir <- "data/processed"
output_dir <- "outputs/macroinvertebrates"
input_all <- file.path(processed_dir, "macro_metrics_v01.rds")
input_study <- file.path(processed_dir, "macro_study_metrics_v01.rds")
for (path in c(input_all, input_study)) {
  if (!file.exists(path)) stop("Missing input: ", path,
                               "\nRun macro_analysis_step01_v02.R first.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(filename) file.path(output_dir, filename)

macro_metrics <- readRDS(input_all)
macro_study_metrics <- readRDS(input_study)
required <- c("site_id", "working_name", "watercourse", "relative_position",
              "site_role", "year", "taxon_richness", "total_abundance")
for (nm in c("macro_metrics", "macro_study_metrics")) {
  missing <- setdiff(required, names(get(nm)))
  if (length(missing)) stop(nm, " missing: ", paste(missing, collapse = ", "))
}
if (!nrow(macro_metrics) || !nrow(macro_study_metrics))
  stop("One or both input datasets contain no samples.")
if (anyNA(macro_study_metrics$year) || anyNA(macro_study_metrics$taxon_richness))
  stop("Missing year or richness in study samples; inspect Script 02 outputs.")
if (any(macro_study_metrics$taxon_richness < 0))
  stop("Negative richness values detected.")

site_order <- c("Upper Gilwiskaw", "B-road Gilwiskaw", "Upper Mease",
                "Restoration Reach Mease", "Birds Hill", "Lower Mease")
long_term_sites <- c("Upper Gilwiskaw", "Restoration Reach Mease", "Lower Mease")
missing_sites <- setdiff(site_order, unique(as.character(macro_study_metrics$working_name)))
if (length(missing_sites))
  stop("Expected study sites missing: ", paste(missing_sites, collapse = ", "))
unexpected_sites <- setdiff(unique(as.character(macro_study_metrics$working_name)),
                           site_order)
if (length(unexpected_sites))
  stop("Unexpected study site labels: ", paste(unexpected_sites, collapse = ", "))

macro_study_metrics <- macro_study_metrics %>%
  mutate(working_name = as.character(working_name))
macro_longterm <- macro_study_metrics %>%
  filter(working_name %in% long_term_sites)
stopifnot(nrow(macro_longterm) > 0)
saveRDS(macro_longterm, file.path(processed_dir, "macro_longterm_v01.rds"))

# ---- 1. Monitoring effort: catchment and study network -------------------
catchment_sampling_effort <- macro_metrics %>%
  count(year, name = "n_samples") %>% arrange(year)
study_sampling_effort <- macro_study_metrics %>%
  count(working_name, year, name = "n_samples") %>%
  arrange(match(working_name, site_order), year)

p_catchment_sampling <- ggplot(catchment_sampling_effort,
                              aes(year, n_samples)) +
  geom_col() +
  labs(title = "Catchment-scale macroinvertebrate sampling effort through time",
       x = "Year", y = "Number of samples") +
  theme_minimal()
ggsave(out("step02_catchment_sampling_effort.png"),
       p_catchment_sampling, width = 8, height = 5, dpi = 300)

p_study_sampling <- study_sampling_effort %>%
  mutate(working_name = factor(working_name, levels = site_order)) %>%
  ggplot(aes(year, n_samples)) +
  geom_col() + facet_wrap(~ working_name, scales = "free_y") +
  labs(title = "Sampling effort at restoration monitoring sites",
       x = "Year", y = "Number of samples") +
  theme_minimal()
ggsave(out("step02_study_sampling_effort.png"),
       p_study_sampling, width = 10, height = 7, dpi = 300)

# ---- 2. Six-site richness summary and presentation table ----------------
study_richness_summary_step02 <- macro_study_metrics %>%
  group_by(working_name, watercourse, relative_position, site_role) %>%
  summarise(
    first_year = min(year), last_year = max(year),
    sample_period = paste0(first_year, "\u2013", last_year),
    n_samples = n(), n_years_sampled = n_distinct(year),
    mean_richness = mean(taxon_richness),
    median_richness = median(taxon_richness),
    min_richness = min(taxon_richness),
    max_richness = max(taxon_richness),
    sd_richness = sd(taxon_richness), .groups = "drop"
  ) %>% arrange(watercourse, relative_position)

study_richness_summary_kable <- macro_study_metrics %>%
  group_by(working_name) %>%
  summarise(
    n = n(), `Sample period` = paste0(min(year), "\u2013", max(year)),
    `Median richness` = median(taxon_richness),
    `Mean richness` = mean(taxon_richness),
    SD = sd(taxon_richness),
    Range = paste0(min(taxon_richness), "\u2013", max(taxon_richness)),
    .groups = "drop"
  ) %>%
  arrange(match(working_name, site_order)) %>%
  select(Site = working_name, everything())
print(kable(study_richness_summary_kable, digits = 1,
            caption = "Taxon richness across the six historical study sites."))

p_study_richness_boxplot <- macro_study_metrics %>%
  mutate(working_name = factor(working_name, levels = rev(site_order))) %>%
  ggplot(aes(working_name, taxon_richness)) +
  geom_boxplot(linewidth = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.35) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = "Study site", y = "Taxon richness per sample") +
  theme_minimal()
ggsave(out("step02_study_richness_boxplot.png"),
       p_study_richness_boxplot, width = 8, height = 5, dpi = 300)

# ---- 3. Long-term site annual summaries and figures ---------------------
longterm_richness_year <- macro_longterm %>%
  group_by(working_name, watercourse, relative_position, site_role, year) %>%
  summarise(
    n_samples = n(),
    mean_richness = mean(taxon_richness),
    median_richness = median(taxon_richness),
    min_richness = min(taxon_richness),
    max_richness = max(taxon_richness),
    .groups = "drop"
  ) %>%
  mutate(working_name = factor(working_name, levels = long_term_sites))

# Preserve the original visual style, including its descriptive LOESS curve.
# Smoothing is illustrative only: uneven sample counts and monitoring years
# prevent interpreting this as a restoration effect or inferential time trend.
p_longterm_richness <- ggplot(longterm_richness_year,
                             aes(year, median_richness)) +
  geom_line(colour = "grey35", linewidth = 0.5) +
  geom_point(aes(size = n_samples), alpha = 0.75) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1.2,
              alpha = 0.15, colour = "#2C6BEA") +
  geom_vline(xintercept = 2021, linetype = "22",
             linewidth = 0.5, colour = "firebrick3") +
  facet_wrap(~ working_name, ncol = 1, scales = "free_y") +
  scale_size_continuous(name = "Annual sample size", range = c(1.5, 5)) +
  scale_x_continuous(breaks = c(1990, 2000, 2010, 2020),
                     expand = expansion(mult = c(0.02, 0.03))) +
  labs(x = "Year", y = "Median taxon richness") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 12),
        panel.spacing = grid::unit(1.2, "lines"),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"))
ggsave(out("step02_longterm_richness_through_time.png"),
       p_longterm_richness, width = 180, height = 230,
       units = "mm", dpi = 300, bg = "white")

p_longterm_richness_boxplot <- ggplot(
  macro_longterm, aes(working_name, taxon_richness)
) +
  geom_boxplot(outlier.alpha = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  coord_flip() +
  labs(title = "Taxon richness at long-term inference sites",
       x = "Site", y = "Taxon richness per sample") +
  theme_minimal()
ggsave(out("step02_longterm_richness_boxplot.png"),
       p_longterm_richness_boxplot, width = 8, height = 5, dpi = 300)

# ---- 4. Abundance diagnostic (not an ecological comparison) -------------
# Exclude zero/negative totals from log10 plotting only; retain all samples
# in the richness analyses and original processed datasets.
diagnostic_data <- macro_study_metrics %>%
  filter(is.finite(total_abundance), total_abundance > 0)
if (nrow(diagnostic_data)) {
  p_richness_abundance <- ggplot(
    diagnostic_data, aes(total_abundance, taxon_richness)
  ) +
    geom_point(alpha = 0.6) +
    scale_x_log10() +
    labs(title = "Diagnostic relationship between richness and recorded abundance",
         subtitle = "Abundance coding systems differ across monitoring years",
         x = "Total recorded abundance (log10 scale)",
         y = "Taxon richness per sample") +
    theme_minimal()
  ggsave(out("step02_richness_abundance_diagnostic.png"),
         p_richness_abundance, width = 8, height = 5, dpi = 300)
} else {
  warning("No positive abundance totals; diagnostic plot not generated.")
}

# ---- 5. Exploratory non-parametric site comparisons ---------------------
# Kruskal-Wallis and BH-adjusted pairwise Wilcoxon tests do not control for
# uneven sampling effort, season, year, or repeated sampling at each site.
richness_kruskal_study <- kruskal.test(
  taxon_richness ~ working_name, data = macro_study_metrics
)
richness_pairwise_study <- pairwise.wilcox.test(
  macro_study_metrics$taxon_richness, macro_study_metrics$working_name,
  p.adjust.method = "BH"
)
richness_kruskal_longterm <- kruskal.test(
  taxon_richness ~ working_name, data = macro_longterm
)
richness_pairwise_longterm <- pairwise.wilcox.test(
  macro_longterm$taxon_richness, macro_longterm$working_name,
  p.adjust.method = "BH"
)

test_summary <- tibble(
  comparison = c("Six study sites", "Three long-term sites"),
  statistic = c(unname(richness_kruskal_study$statistic),
                unname(richness_kruskal_longterm$statistic)),
  df = c(unname(richness_kruskal_study$parameter),
         unname(richness_kruskal_longterm$parameter)),
  p_value = c(richness_kruskal_study$p.value,
              richness_kruskal_longterm$p.value)
)
pairwise_to_csv <- function(test, label) {
  as.data.frame(as.table(test$p.value), stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    rename(site_1 = Var1, site_2 = Var2, p_adjusted_BH = Freq) %>%
    filter(!is.na(p_adjusted_BH)) %>%
    mutate(comparison = label, .before = 1)
}
pairwise_summary <- bind_rows(
  pairwise_to_csv(richness_pairwise_study, "Six study sites"),
  pairwise_to_csv(richness_pairwise_longterm, "Three long-term sites")
)

# ---- 6. Descriptive monitoring-period context ---------------------------
# 2021 is the implementation-year boundary, not a verified completion date
# for every individual sample; do not interpret as a causal intervention test.
macro_longterm_period <- macro_longterm %>%
  mutate(restoration_period = case_when(
    year < 2021 ~ "Historic baseline",
    year >= 2021 ~ "Contemporary monitoring period",
    TRUE ~ NA_character_
  ))
period_richness_summary <- macro_longterm_period %>%
  group_by(working_name, restoration_period) %>%
  summarise(
    n_samples = n(), mean_richness = mean(taxon_richness),
    median_richness = median(taxon_richness),
    min_richness = min(taxon_richness),
    max_richness = max(taxon_richness),
    .groups = "drop"
  )
p_period_richness <- ggplot(
  macro_longterm_period, aes(restoration_period, taxon_richness)
) +
  geom_boxplot(outlier.alpha = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  facet_wrap(~ working_name) +
  labs(title = "Taxon richness by monitoring period",
       subtitle = "Descriptive comparison only; sampling is uneven across periods",
       x = NULL, y = "Taxon richness per sample") +
  theme_minimal()
ggsave(out("step02_pre_post_richness_context.png"),
       p_period_richness, width = 10, height = 6, dpi = 300)

macro_longterm_year_summary <- macro_longterm %>%
  group_by(year) %>%
  summarise(n_samples = n(),
            mean_richness = mean(taxon_richness),
            median_richness = median(taxon_richness),
            .groups = "drop")

# ---- 7. Export all summaries --------------------------------------------
write_csv(catchment_sampling_effort, out("step02_catchment_sampling_effort.csv"))
write_csv(study_sampling_effort, out("step02_study_sampling_effort.csv"))
write_csv(study_richness_summary_step02, out("step02_study_richness_summary.csv"))
write_csv(study_richness_summary_kable, out("step02_study_richness_table.csv"))
write_csv(longterm_richness_year, out("step02_longterm_richness_year.csv"))
write_csv(period_richness_summary, out("step02_pre_post_richness_summary.csv"))
write_csv(test_summary, out("step02_richness_kruskal_tests.csv"))
write_csv(pairwise_summary, out("step02_richness_pairwise_BH.csv"))
write_csv(macro_longterm_year_summary, out("step02_longterm_year_summary.csv"))

# ---- 8. End-of-script checks --------------------------------------------
stopifnot(sum(catchment_sampling_effort$n_samples) == nrow(macro_metrics))
stopifnot(sum(study_sampling_effort$n_samples) == nrow(macro_study_metrics))
stopifnot(sum(study_richness_summary_step02$n_samples) ==
            nrow(macro_study_metrics))
message("Script 03 complete: macro_analysis_step02_v02.R")
message("Catchment samples: ", nrow(macro_metrics),
        "; study samples: ", nrow(macro_study_metrics),
        "; long-term samples: ", nrow(macro_longterm))
message("Six-site Kruskal-Wallis p: ",
        signif(richness_kruskal_study$p.value, 4))
message("Long-term Kruskal-Wallis p: ",
        signif(richness_kruskal_longterm$p.value, 4))
message("Outputs: ", output_dir)
