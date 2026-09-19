# ============================================================================
# River Mease floodplain restoration | Macroinvertebrates
# Script 02: Historical taxon richness and monitoring coverage
# File: macro_analysis_step01_v02.R
# ============================================================================
# PURPOSE
#   Calculate abundance-independent taxon richness for every historical sample,
#   identify the six study sites and three long-term sites, and export monitoring
#   coverage, descriptive richness summaries and exploratory figures.
#
# RUN ORDER / WORKING DIRECTORY
#   Run AFTER macro_historical_v03.R, from the EA_data project root (not scripts/).
#   Example: setwd('/home/hamish/Documents/NTU/research_project/data_analysis/data/EA_data')
#
# INPUTS (created by script 01)
#   data/processed/macro_clean_v01.rds
#   data/processed/macro_taxa_cols_v01.rds
#
# OUTPUTS (v01 filenames retained for downstream compatibility)
#   data/processed/macro_metrics_v01.rds
#   data/processed/macro_study_metrics_v01.rds
#   data/processed/macro_longterm_v01.rds
#   outputs/macroinvertebrates/study_inventory.csv
#   outputs/macroinvertebrates/study_richness_summary.csv
#   outputs/macroinvertebrates/catchment_richness_through_time.png
#   outputs/macroinvertebrates/study_site_richness_through_time.png
#   Additional monitoring-coverage CSVs and figures documented below.
#
# ANALYTICAL DECISIONS
#   Historical abundance recording systems differ. Richness counts taxa with
#   recorded values > 0; missing values do not count as detections. Summed
#   abundance is saved ONLY as a diagnostic, not compared across years.
#   Taxon richness reflects the source taxonomic columns, not a harmonised
#   BMWP-family metric. Later scripts calculate BMWP/ASPT separately.
#   Site labels are descriptive; historic and contemporary sampling locations
#   should not be assumed to be geographically identical.
#
# STATUS: Revised from macro_analysis_step01_v01.R; requires execution and
# comparison with the submitted analyses before publication as validated code.
# ============================================================================

# ---- 0. Packages, paths and input checks ---------------------------------
library(tidyverse)
library(knitr)

processed_dir <- 'data/processed'
output_dir <- 'outputs/macroinvertebrates'
input_macro <- file.path(processed_dir, 'macro_clean_v01.rds')
input_taxa <- file.path(processed_dir, 'macro_taxa_cols_v01.rds')

for (path in c(input_macro, input_taxa)) {
  if (!file.exists(path)) stop('Missing input: ', path,
                               '\nRun macro_historical_v03.R from the project root first.')
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Read and validate the cleaned dataset ----------------------------
macro <- readRDS(input_macro)
taxa_cols <- readRDS(input_taxa)
required_cols <- c('site_id', 'sample_id', 'year', 'working_name',
                   'watercourse', 'relative_position', 'site_role')
missing_cols <- setdiff(required_cols, names(macro))
if (length(missing_cols)) stop('Missing metadata columns: ',
                               paste(missing_cols, collapse = ', '))
missing_taxa <- setdiff(taxa_cols, names(macro))
if (length(missing_taxa)) stop('Missing taxa columns: ',
                               paste(missing_taxa, collapse = ', '))
if (!length(taxa_cols)) stop('Taxa column list is empty.')
if (anyDuplicated(taxa_cols)) stop('Duplicate names in taxa_cols.')
if (!all(vapply(macro[taxa_cols], is.numeric, logical(1)))) {
  stop('Non-numeric taxa columns found; check the historical import.')
}
if (anyNA(macro$site_id)) warning('Missing site IDs detected.')

# ---- 2. Define the six-site study network --------------------------------
# Names and role descriptions match the original Step 01 analysis. Script 01
# already adds similarly named columns: do NOT join another register without
# removing/reconciling those columns, as this creates .x/.y suffixes.
study_site_lookup <- tribble(
  ~site_id, ~study_name, ~study_watercourse, ~study_position, ~study_role,
  '50427', 'Upper Gilwiskaw', 'Gilwiskaw', 'Upstream', 'Long-term upstream reference',
  '211791', 'B-road Gilwiskaw', 'Gilwiskaw', 'Upstream', 'Contemporary upstream context',
  '156912', 'Upper Mease', 'Mease', 'Upstream', 'Supplementary upstream context',
  '52073', 'Restoration Reach Mease', 'Mease', 'Restoration', 'Long-term restoration site',
  '50447', 'Lower Mease', 'Mease', 'Downstream', 'Long-term downstream response',
  'Birds Hill', 'Birds Hill', 'Mease', 'Downstream', 'Contemporary downstream context'
)
stopifnot(!anyDuplicated(study_site_lookup$site_id))

# ---- 3. Calculate per-sample metrics -------------------------------------
# Matrix operations are equivalent to the original rowwise sums, but avoid
# creating hundreds of thousands of intermediate rowwise operations.
taxa_matrix <- as.matrix(macro[, taxa_cols, drop = FALSE])
if (any(taxa_matrix < 0, na.rm = TRUE)) {
  warning('Negative taxa values detected; inspect source before interpreting richness.')
}
macro_metrics <- macro %>%
  mutate(
    site_id = as.character(site_id),
    taxon_richness = rowSums(taxa_matrix > 0, na.rm = TRUE),
    total_abundance = rowSums(taxa_matrix, na.rm = TRUE),
    n_taxa_missing = rowSums(is.na(taxa_matrix))
  )
# total_abundance is a diagnostic only: historical counting methods differ.

# ---- 4. Select study and long-term samples -------------------------------
# Replace Script 01's general-purpose site metadata with the Step 01 analysis
# register explicitly. Avoid accidental duplicate working_name columns.
macro_study_metrics <- macro_metrics %>%
  filter(site_id %in% study_site_lookup$site_id) %>%
  select(-any_of(c('working_name', 'watercourse', 'relative_position', 'site_role'))) %>%
  left_join(study_site_lookup, by = 'site_id') %>%
  rename(working_name = study_name,
         watercourse = study_watercourse,
         relative_position = study_position,
         site_role = study_role)

long_term_sites <- c('Upper Gilwiskaw', 'Restoration Reach Mease', 'Lower Mease')
macro_longterm <- macro_study_metrics %>%
  filter(working_name %in% long_term_sites)

missing_sites <- setdiff(study_site_lookup$site_id, unique(macro_study_metrics$site_id))
if (length(missing_sites)) warning('Study site IDs absent from data: ',
                                   paste(missing_sites, collapse = ', '))
if (nrow(macro_study_metrics) != sum(macro_metrics$site_id %in% study_site_lookup$site_id)) {
  stop('Unexpected row duplication or loss after joining the study-site register.')
}

# ---- 5. Monitoring coverage ---------------------------------------------
study_inventory <- macro_study_metrics %>%
  group_by(site_id, working_name, watercourse, relative_position, site_role) %>%
  summarise(first_year = min(year, na.rm = TRUE),
            last_year = max(year, na.rm = TRUE),
            n_samples = n(),
            n_years_sampled = n_distinct(year, na.rm = TRUE),
            .groups = 'drop') %>%
  arrange(watercourse, first_year) %>%
  transmute(Site = working_name, Role = site_role,
            `First year` = first_year, `Last year` = last_year,
            Samples = n_samples, `Monitoring years` = n_years_sampled)
print(knitr::kable(study_inventory,
                   caption = 'Macroinvertebrate monitoring coverage at the six study sites.'))
write_csv(study_inventory, file.path(output_dir, 'study_inventory.csv'))

catchment_sampling_effort <- macro_metrics %>%
  count(year, name = 'n_samples') %>% arrange(year)
study_sampling_effort <- macro_study_metrics %>%
  count(working_name, year, name = 'n_samples') %>% arrange(working_name, year)
write_csv(catchment_sampling_effort,
          file.path(output_dir, 'step01_catchment_sampling_effort.csv'))
write_csv(study_sampling_effort,
          file.path(output_dir, 'step01_study_sampling_effort.csv'))

p_catchment_sampling <- ggplot(catchment_sampling_effort, aes(year, n_samples)) +
  geom_col() + labs(title = 'Macroinvertebrate sampling effort through time',
                    x = 'Year', y = 'Number of samples') + theme_minimal()
ggsave(file.path(output_dir, 'step01_catchment_sampling_effort.png'),
       p_catchment_sampling, width = 8, height = 5, dpi = 300)

p_study_sampling <- ggplot(study_sampling_effort, aes(year, n_samples)) +
  geom_col() + facet_wrap(~working_name, scales = 'free_y') +
  labs(title = 'Sampling effort at restoration monitoring sites',
       x = 'Year', y = 'Number of samples') + theme_minimal()
ggsave(file.path(output_dir, 'step01_study_sampling_effort.png'),
       p_study_sampling, width = 10, height = 7, dpi = 300)

# ---- 6. Catchment richness through time ---------------------------------
catchment_richness_year <- macro_metrics %>%
  group_by(year) %>%
  summarise(n_samples = n(),
            mean_richness = mean(taxon_richness, na.rm = TRUE),
            median_richness = median(taxon_richness, na.rm = TRUE),
            .groups = 'drop') %>% arrange(year)
write_csv(catchment_richness_year,
          file.path(output_dir, 'step01_catchment_richness_year.csv'))

p_catchment_richness <- ggplot(catchment_richness_year,
                              aes(year, median_richness)) +
  geom_line() + geom_point() +
  labs(title = 'Catchment-scale macroinvertebrate richness through time',
       x = 'Year', y = 'Median taxon richness per sample') + theme_minimal()
ggsave(file.path(output_dir, 'catchment_richness_through_time.png'),
       p_catchment_richness, width = 8, height = 5, dpi = 300)

# ---- 7. Study-site richness through time ---------------------------------
# LOESS is retained from the original exploratory script; it is descriptive
# only and should not be interpreted as a restoration effect or formal trend.
p_study_richness <- ggplot(macro_study_metrics,
                          aes(year, taxon_richness)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = 'loess', se = TRUE) +
  facet_wrap(~working_name, scales = 'free_y') +
  labs(title = 'Taxon richness through time at restoration-network sites',
       x = 'Year', y = 'Taxon richness per sample') + theme_minimal()
ggsave(file.path(output_dir, 'study_site_richness_through_time.png'),
       p_study_richness, width = 10, height = 7, dpi = 300)

# ---- 8. Study-site richness summary -------------------------------------
study_richness_summary <- macro_study_metrics %>%
  group_by(working_name, site_id, watercourse, relative_position, site_role) %>%
  summarise(first_year = min(year, na.rm = TRUE),
            last_year = max(year, na.rm = TRUE),
            n_samples = n(),
            n_years_sampled = n_distinct(year, na.rm = TRUE),
            mean_richness = mean(taxon_richness, na.rm = TRUE),
            median_richness = median(taxon_richness, na.rm = TRUE),
            min_richness = min(taxon_richness, na.rm = TRUE),
            max_richness = max(taxon_richness, na.rm = TRUE),
            mean_total_abundance = mean(total_abundance, na.rm = TRUE),
            median_total_abundance = median(total_abundance, na.rm = TRUE),
            .groups = 'drop') %>%
  arrange(watercourse, relative_position)
# Abundance columns above are diagnostic only; do not compare between years.
write_csv(study_richness_summary,
          file.path(output_dir, 'study_richness_summary.csv'))

# ---- 9. Save objects needed by downstream scripts -----------------------
saveRDS(macro_metrics, file.path(processed_dir, 'macro_metrics_v01.rds'))
saveRDS(macro_study_metrics,
        file.path(processed_dir, 'macro_study_metrics_v01.rds'))
saveRDS(macro_longterm, file.path(processed_dir, 'macro_longterm_v01.rds'))

# ---- 10. End-of-script checks -------------------------------------------
stopifnot(nrow(macro_metrics) == nrow(macro),
          nrow(macro_longterm) <= nrow(macro_study_metrics),
          all(macro_study_metrics$site_id %in% study_site_lookup$site_id),
          all(macro_metrics$taxon_richness >= 0))
message('Step 02 complete: ', nrow(macro_metrics), ' catchment samples; ',
        nrow(macro_study_metrics), ' study-site samples; ',
        nrow(macro_longterm), ' long-term-site samples.')
message('Outputs saved to: ', processed_dir, ' and ', output_dir)
