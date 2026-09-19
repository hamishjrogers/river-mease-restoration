# ============================================================================
# River Mease | Macroinvertebrates | Script 07
# macro_analysis_step06_v02.R — combine historic and 2026 BMWP/ASPT
# ============================================================================
# Run from EA_data project root, after macro_analysis_step04_v05.R.
# Historical scoring is retained unchanged; six contemporary samples are
# scored from their original operational-group field matrix. Raw taxon richness
# is NOT compared between historic and contemporary samples because taxonomic
# resolution differs. Two 2026 locations are approximate historic-site matches,
# not identical sampling coordinates. No causal restoration inference is made.
# Inputs: data/processed/macro_bmwp_scored_v01.rds;
#         data/2026_measemacrosampling.csv
# Outputs: original combined BMWP and 2026 community RDS/CSV filenames.
# The date below is retained from the original script; confirm against field log.
# ============================================================================
library(tidyverse)

processed_dir <- 'data/processed'
input_historic <- file.path(processed_dir, 'macro_bmwp_scored_v01.rds')
input_2026 <- 'data/2026_measemacrosampling.csv'
for (p in c(input_historic, input_2026)) {
  if (!file.exists(p)) stop('Missing input: ', p)
}
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
survey_date_2026 <- as.Date('2026-07-19') # Original script date; verify with field log.
site_columns_2026 <- c('Upstream Gilwiskaw', 'Wetland Gilwiskaw',
                       'Upstream Mease', 'Wetland Mease', 'Confluence',
                       'Downstream Mease')
required_2026 <- c('Taxa', 'BMWP Score', site_columns_2026)

# ---- 1. Read and validate historical and contemporary inputs ------------
macro_bmwp_historic <- readRDS(input_historic)
required_historic <- c('working_name', 'watercourse', 'relative_position',
                       'site_role', 'sample_date', 'year', 'taxon_richness',
                       'n_bmwp_groups', 'bmwp', 'aspt')
missing_historic <- setdiff(required_historic, names(macro_bmwp_historic))
if (length(missing_historic)) stop('Historical input missing: ',
                                  paste(missing_historic, collapse = ', '))
if (nrow(macro_bmwp_historic) != 147L ||
    n_distinct(macro_bmwp_historic$site_id) != 6L)
  stop('Historical input does not contain the expected 147 samples / 6 sites.')
if (anyNA(macro_bmwp_historic$bmwp) || anyNA(macro_bmwp_historic$aspt))
  stop('Historical BMWP/ASPT contains missing values; inspect scoring output.')

bmwp_2026_raw <- read_csv(input_2026, show_col_types = FALSE)
missing_2026 <- setdiff(required_2026, names(bmwp_2026_raw))
if (length(missing_2026)) stop('2026 field sheet missing: ',
                               paste(missing_2026, collapse = ', '))

# ---- 2. Identify scoring rows and preserve contemporary occurrences -----
# Section headings/metadata have no numeric BMWP score and are excluded.
# Abundance parsing matches the original script: blank/NA/N/A/- -> zero.
score_text <- trimws(as.character(bmwp_2026_raw[['BMWP Score']]))
score_num <- suppressWarnings(parse_double(score_text,
                                            na = c('', 'NA', 'N/A', '-')))
scoring_rows <- !is.na(score_num)
if (!any(scoring_rows)) stop('No numeric BMWP scoring rows found.')
if (any(is.na(bmwp_2026_raw$Taxa[scoring_rows]) |
        !nzchar(trimws(as.character(bmwp_2026_raw$Taxa[scoring_rows])))))
  stop('One or more scoring rows have a blank operational-group name.')
if (any(score_num[scoring_rows] < 0)) stop('Negative BMWP scores detected.')

bmwp_2026_long <- bmwp_2026_raw %>%
  mutate(bmwp_score_numeric = score_num) %>%
  filter(!is.na(bmwp_score_numeric)) %>%
  pivot_longer(cols = all_of(site_columns_2026),
               names_to = 'field_site', values_to = 'abundance_raw') %>%
  mutate(abundance = suppressWarnings(parse_double(
           trimws(as.character(abundance_raw)), na = c('', 'NA', 'N/A', '-'))),
         abundance = replace_na(abundance, 0),
         present = abundance > 0,
         `BMWP Score` = bmwp_score_numeric) %>%
  select(-bmwp_score_numeric)

# A group must have one row/score in the source matrix. A repeated group
# would inflate the score even if its rows carry the same score.
duplicate_groups <- bmwp_2026_long %>%
  count(field_site, Taxa, name = 'n') %>% filter(n > 1)
if (nrow(duplicate_groups)) {
  print(duplicate_groups, n = Inf)
  stop('Repeated operational-group rows in 2026 scoring sheet.')
}
if (anyNA(bmwp_2026_long$abundance)) stop('Unresolved 2026 abundance values.')
if (any(bmwp_2026_long$abundance < 0)) stop('Negative abundance detected.')

# ---- 3. Calculate contemporary BMWP/ASPT -------------------------------
bmwp_2026_scores <- bmwp_2026_long %>%
  filter(present) %>%
  distinct(field_site, Taxa, `BMWP Score`, .keep_all = TRUE) %>%
  group_by(field_site) %>%
  summarise(n_bmwp_groups = n_distinct(Taxa),
            bmwp = sum(`BMWP Score`),
            aspt = bmwp / n_bmwp_groups, .groups = 'drop')
# Retain every site even if it has zero scoring groups.
bmwp_2026_scores <- tibble(field_site = site_columns_2026) %>%
  left_join(bmwp_2026_scores, by = 'field_site') %>%
  mutate(n_bmwp_groups = replace_na(n_bmwp_groups, 0L),
         bmwp = replace_na(bmwp, 0))
if (nrow(bmwp_2026_scores) != 6L ||
    anyNA(bmwp_2026_scores$aspt))
  stop('Expected six 2026 sites with defined ASPT; check source scoring.')

# ---- 4. Harmonise site metadata (original mapping retained) ------------
bmwp_2026_harmonised <- bmwp_2026_scores %>%
  mutate(
    working_name = case_when(
      field_site == 'Upstream Gilwiskaw' ~ 'B-road Gilwiskaw',
      field_site == 'Wetland Mease' ~ 'Restoration Reach Mease',
      TRUE ~ field_site),
    watercourse = case_when(
      field_site %in% c('Upstream Gilwiskaw', 'Wetland Gilwiskaw') ~ 'Gilwiskaw',
      field_site %in% c('Upstream Mease', 'Wetland Mease', 'Downstream Mease') ~ 'Mease',
      field_site == 'Confluence' ~ 'Mease-Gilwiskaw confluence'),
    relative_position = case_when(
      field_site %in% c('Upstream Gilwiskaw', 'Upstream Mease') ~ 'Upstream',
      field_site %in% c('Wetland Gilwiskaw', 'Wetland Mease') ~ 'Restoration',
      field_site == 'Confluence' ~ 'Confluence',
      field_site == 'Downstream Mease' ~ 'Downstream'),
    site_role = case_when(
      field_site == 'Upstream Gilwiskaw' ~ 'Contemporary upstream context',
      field_site == 'Wetland Gilwiskaw' ~ 'Contemporary Gilwiskaw wetland site',
      field_site == 'Upstream Mease' ~ 'Contemporary Mease upstream site',
      field_site == 'Wetland Mease' ~ 'Long-term restoration site',
      field_site == 'Confluence' ~ 'Contemporary confluence site',
      field_site == 'Downstream Mease' ~ 'Contemporary downstream site'),
    date_sampled = survey_date_2026,
    year = 2026,
    data_source = '2026 field survey')

historic_names <- unique(as.character(macro_bmwp_historic$working_name))
actual_matches <- bmwp_2026_harmonised %>%
  filter(working_name %in% historic_names) %>% pull(working_name)
expected_matches <- c('B-road Gilwiskaw', 'Restoration Reach Mease')
if (!setequal(actual_matches, expected_matches))
  stop('Unexpected historic/2026 site-name matches: ',
       paste(actual_matches, collapse = ', '))

# ---- 5. Create deliberately compact, comparable metrics table ----------
bmwp_historic_for_join <- macro_bmwp_historic %>%
  transmute(field_site = NA_character_,
            working_name = as.character(working_name),
            watercourse = as.character(watercourse),
            relative_position = as.character(relative_position),
            site_role = as.character(site_role),
            sample_date = as.Date(sample_date),
            year = as.numeric(year),
            taxon_richness = as.integer(taxon_richness),
            n_bmwp_groups = as.integer(n_bmwp_groups),
            bmwp = as.numeric(bmwp), aspt = as.numeric(aspt),
            data_source = 'Historic secondary data')
bmwp_2026_for_join <- bmwp_2026_harmonised %>%
  transmute(field_site, working_name, watercourse, relative_position,
            site_role, sample_date = as.Date(date_sampled),
            year = as.numeric(year), taxon_richness = NA_integer_,
            n_bmwp_groups = as.integer(n_bmwp_groups),
            bmwp = as.numeric(bmwp), aspt = as.numeric(aspt), data_source)
macro_bmwp_combined <- bind_rows(bmwp_historic_for_join,
                                bmwp_2026_for_join) %>%
  arrange(working_name, year, sample_date)
stopifnot(nrow(macro_bmwp_combined) == 153L,
          sum(macro_bmwp_combined$year == 2026) == 6L,
          !anyNA(macro_bmwp_combined$bmwp),
          !anyNA(macro_bmwp_combined$aspt),
          all(is.na(bmwp_2026_for_join$taxon_richness)))

# ---- 6. Export: original filenames preserved ----------------------------
combined_rds_path <- file.path(processed_dir,
                               'macro_bmwp_scored_combined_2026_v01.rds')
combined_csv_path <- file.path(processed_dir,
                               'macro_bmwp_scored_combined_2026_v01.csv')
saveRDS(macro_bmwp_combined, combined_rds_path)
write_csv(macro_bmwp_combined, combined_csv_path)
saveRDS(bmwp_2026_long,
        file.path(processed_dir, 'bmwp_2026_community_long_v01.rds'))
write_csv(bmwp_2026_long,
          file.path(processed_dir, 'bmwp_2026_community_long_v01.csv'))

# ---- 7. Console audit ---------------------------------------------------
cat('\n=== Historic + 2026 BMWP/ASPT combination complete ===\n')
print(macro_bmwp_combined %>% count(data_source), n = Inf)
print(bmwp_2026_harmonised %>%
        select(field_site, working_name, n_bmwp_groups, bmwp, aspt),
      n = Inf, width = Inf)
cat('Combined rows:', nrow(macro_bmwp_combined),
    '| Historical:', nrow(macro_bmwp_historic),
    '| Contemporary:', nrow(bmwp_2026_harmonised), '\n')
cat('2026 survey date retained from original script:',
    as.character(survey_date_2026), '(verify against field log)\n')
cat('Saved:', combined_rds_path, '\n')
