# ---- Script information --------------------------------------------------
# Title: macro_historical_v03.R
# Author: HR (restructured from macro_historic_v01.R)
# Purpose: Import and quality-check the transposed historic macroinvertebrate
# workbook; identify metadata/taxa fields; describe monitoring coverage;
# attach the six-site study register; export cleaned data and site inventories.
#
# Run order: 1 (before macro_analysis_step01_v01.R and all subsequent analyses).
# Run from: River Mease project root (relative paths are used throughout).
# Input: data/Abundance Database Final.xlsx, worksheet "Final Main".
# Outputs: data/processed/macro_clean_v01.rds
#          data/processed/macro_core_historical_v01.rds
#          data/processed/macro_study_sites_v01.rds
#          data/processed/macro_taxa_cols_v01.rds
#          outputs/macroinvertebrates/*.csv (inventories and site register)
#
# Analytical decisions and limitations:
# - Workbook orientation is unusual: source variable names are in column 1,
#   and each remaining source column represents one sample; transpose once.
# - Taxa columns retain the original recorded abundance values. Do not compare
#   abundance across years: source recording systems changed over time.
#   Downstream long-term analyses use presence/absence where appropriate.
# - Missing taxa entries remain NA in the cleaned source data; downstream
#   analyses must document any absence/zero convention they apply.
# - The six-site register identifies the study network; all other source sites
#   are retained in the full catchment dataset as context.
# - A study-site label does not establish exact geographic equivalence between
#   historic and 2026 field locations; interpret comparisons accordingly.
# - This revision preserves the original v01 output filenames for compatibility
#   with downstream scripts. It does not silently alter historical results.
#
# Reproducibility: input workbook is not distributed in the public repository.
# Status: code reviewed; requires execution against the source workbook and
# downstream output comparison before being described as validated.

# ---- 0. Packages, input and output paths --------------------------------

library(tidyverse)
library(readxl)
library(janitor)
library(lubridate)
library(rnrfa) # osg_parse() for optional QGIS coordinates

setwd("/home/hamish/Documents/NTU/research_project/data_analysis/data/EA_data")
input_file <- "data/Abundance Database Final.xlsx"
output_dir <- "outputs/macroinvertebrates"
processed_dir <- "data/processed"
if (!file.exists(input_file)) stop("Missing source workbook: ", input_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)


# ---- 1. Import raw workbook --------------------------------------------

# Read without column headings because the source sheet stores field names
# vertically in its first column.
macro_raw <- read_excel(
  input_file,
  sheet = "Final Main",
  col_names = FALSE
)


# ---- 2. Transpose dataset ----------------------------------------------

# First column contains the variable names
variable_names <- macro_raw[[1]]

# Replace any blank or missing variable names
variable_names <- ifelse(
  is.na(variable_names) | variable_names == "",
  paste0("unknown_", seq_along(variable_names)),
  variable_names
)

# Remove first column, transpose, convert to tibble
macro <- macro_raw %>%
  select(-1) %>%
  t() %>%
  as_tibble(.name_repair = "minimal")

# Apply variable names
names(macro) <- variable_names

# Clean the column names
names(macro) <- make_clean_names(names(macro))
if (anyDuplicated(names(macro))) stop("Duplicate column names after cleaning")
required_metadata <- c("site_id", "sample_id", "sample_location",
                       "sample_ngr", "sample_date", "data_source", "method")
missing_metadata <- setdiff(required_metadata, names(macro))
if (length(missing_metadata)) stop("Missing metadata columns: ",
                                  paste(missing_metadata, collapse = ", "))

# Check
names(macro)[grepl("unknown", names(macro))]
# [1] "unknown_550" "unknown_551" "unknown_552" "unknown_553" "unknown_554" "unknown_555" "unknown_556" "unknown_557"

# Drop blank/unknown rows from the original sheet
macro <- macro %>%
  select(-starts_with("unknown_"))

# ---- 3. Basic variable formatting --------------------------------------

macro <- macro %>%
  mutate(
    sample_date = as.Date(as.numeric(sample_date), origin = "1899-12-30"),
    year = year(sample_date),
    data_source = as.factor(data_source),
    site_id = as.factor(site_id),
    sample_location = as.factor(sample_location),
    method = as.factor(method)
  )


# ---- 4. Identify metadata and taxa columns -----------------------------

metadata_cols <- c(
  "site_id",
  "sample_id",
  "sample_location",
  "data_source",
  "sample_ngr",
  "sample_date",
  "method",
  "unit",
  "year",
  "working_name",
  "watercourse",
  "relative_position",
  "wq_sample_location",
  "site_role"
)

taxa_cols <- setdiff(names(macro), metadata_cols)
if (!length(taxa_cols)) stop("No taxa columns found after metadata separation")


# Convert taxa columns to numeric
macro <- macro %>%
  mutate(across(all_of(taxa_cols), ~ suppressWarnings(as.numeric(.x))))

# QC: report source entries that are nonblank but could not be converted.
# Do not reinterpret these entries as absences without checking the workbook.
non_numeric_taxa <- macro_raw %>%
  select(-1) %>%
  t() %>%
  as_tibble(.name_repair = "minimal")
names(non_numeric_taxa) <- make_clean_names(variable_names)
non_numeric_taxa <- non_numeric_taxa %>% select(all_of(taxa_cols))
conversion_issues <- tibble(
  taxon = taxa_cols,
  n_non_numeric = vapply(taxa_cols, function(nm) {
    original <- trimws(as.character(non_numeric_taxa[[nm]]))
    sum(!is.na(original) & nzchar(original) &
          is.na(suppressWarnings(as.numeric(original))))
  }, integer(1))
) %>% filter(n_non_numeric > 0)
if (nrow(conversion_issues)) {
  warning("Non-numeric taxa entries detected; inspect taxa_conversion_issues.csv")
  write_csv(conversion_issues, file.path(output_dir, "taxa_conversion_issues.csv"))
}


# ---- 5. Dataset inventory ----------------------------------------------

dataset_inventory <- tibble(
  n_samples = nrow(macro),
  n_sites = n_distinct(macro$site_id),
  n_taxa_columns = length(taxa_cols),
  first_year = min(macro$year, na.rm = TRUE),
  last_year = max(macro$year, na.rm = TRUE),
  data_sources = paste(unique(macro$data_source), collapse = ", ")
)

dataset_inventory


# ---- 6. Samples by year -------------------------------------------------

samples_by_year <- macro %>%
  count(year, name = "n_samples") %>%
  arrange(year)

samples_by_year


# ---- 7. Samples by site -------------------------------------------------

samples_by_site <- macro %>%
  count(site_id, name = "n_samples") %>%
  arrange(desc(n_samples))

samples_by_site


# ---- 8. Samples by site and location -----------------------------------

samples_by_site_location <- macro %>%
  count(site_id, sample_location, name = "n_samples") %>%
  arrange(site_id, sample_location)

samples_by_site_location


# ---- 9. Site summary ----------------------------------------------------

site_summary <- macro %>%
  group_by(site_id, sample_location, sample_ngr) %>%
  summarise(
    first_year = min(year, na.rm = TRUE),
    last_year = max(year, na.rm = TRUE),
    n_samples = n(),
    n_years_sampled = n_distinct(year),
    data_sources = paste(unique(data_source), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(site_id, sample_location)

site_summary

# ---- 9b. Export site summary for QGIS ----------------------------------

site_summary_qgis <- site_summary %>%
  mutate(
    sample_ngr = str_replace_all(sample_ngr, "\\s+", "")
  ) %>%
  
  rowwise() %>%
  mutate(
    coords = list(tryCatch(osg_parse(sample_ngr),
                           error = function(e) list(easting = NA_real_,
                                                    northing = NA_real_))),
    easting = coords$easting,
    northing = coords$northing
  ) %>%
  ungroup() %>%
  
  select(
    site_id,
    sample_location,
    sample_ngr,
    easting,
    northing,
    first_year,
    last_year,
    n_samples,
    n_years_sampled,
    data_sources
  ) %>%
  
  arrange(site_id, sample_location)

write_csv(
  site_summary_qgis,
  "outputs/macroinvertebrates/site_summary_qgis.csv"
)

# ---- 9c. Master site register ------------------------------------------

site_register <- tribble(
  ~site_id,       ~working_name,          ~watercourse,        ~relative_position,                  ~wq_sample_location, ~site_role,
  "50427",        "Upper Gilwiskaw",      "Gilwiskaw Brook",   "Upstream Gilwiskaw - furthest",     "No",                "Core historical",
  "211791",       "B-road Gilwiskaw",     "Gilwiskaw Brook",   "Upstream Gilwiskaw - B road",       "Yes",               "Contemporary integration",
  "52073",        "Restoration Reach",    "River Mease",       "Restoration reach - Mease",         "Yes",               "Core historical",
  "Birds Hill",   "Birds Hill",           "River Mease",       "Downstream Mease",                  "No",                "Contemporary integration",
  "50447",        "Lower Mease",          "River Mease",       "Downstream Mease - further",        "No",                "Core historical",
  "156912",       "Upstream Mease",        "River Mease",       "Upstream Mease",                    "No",                "Supplementary context"
)

# ---- 9d. Add site roles to macro dataset --------------------------------

macro <- macro %>%
  mutate(site_id = as.character(site_id)) %>%
  left_join(site_register, by = "site_id") %>%
  mutate(
    site_role = replace_na(site_role, "Catchment context"),
    working_name = coalesce(working_name, site_id),
    watercourse = replace_na(watercourse, "Unclassified"),
    relative_position = replace_na(relative_position, "Wider catchment"),
    wq_sample_location = replace_na(wq_sample_location, "No")
  )

# ---- 10. Site-year summary ---------------------------------------------

site_year_summary <- macro %>%
  count(site_id, sample_location, year, name = "n_samples") %>%
  arrange(site_id, sample_location, year)

site_year_summary


# ---- 11. Check abundance value types -----------------------------------

abundance_values_by_year <- macro %>%
  select(year, all_of(taxa_cols)) %>%
  pivot_longer(
    cols = all_of(taxa_cols),
    names_to = "taxon",
    values_to = "abundance"
  ) %>%
  filter(!is.na(abundance), abundance > 0) %>%
  distinct(year, abundance) %>%
  arrange(year, abundance)

abundance_values_by_year


# ---- 12. Export outputs -------------------------------------------------

write_csv(dataset_inventory, "outputs/macroinvertebrates/dataset_inventory.csv")
write_csv(samples_by_year, "outputs/macroinvertebrates/samples_by_year.csv")
write_csv(samples_by_site, "outputs/macroinvertebrates/samples_by_site.csv")
write_csv(samples_by_site_location, "outputs/macroinvertebrates/samples_by_site_location.csv")
write_csv(site_summary, "outputs/macroinvertebrates/site_summary.csv")
write_csv(site_year_summary, "outputs/macroinvertebrates/site_year_summary.csv")
write_csv(abundance_values_by_year, "outputs/macroinvertebrates/abundance_values_by_year.csv")
write_csv(
  site_register,
  "outputs/macroinvertebrates/site_register.csv"
)

macro_core <- macro %>%
  filter(site_role == "Core historical")

macro_study <- macro %>%
  filter(site_role %in% c(
    "Core historical",
    "Contemporary integration",
    "Supplementary context"
  ))

saveRDS(macro, "data/processed/macro_clean_v01.rds")
saveRDS(macro_core, "data/processed/macro_core_historical_v01.rds")
saveRDS(macro_study, "data/processed/macro_study_sites_v01.rds")
saveRDS(taxa_cols, "data/processed/macro_taxa_cols_v01.rds")


# ---- 13. End-of-script quality-control summary --------------------------
# Print the quantities used to verify this import against the original run.
# Differences require investigation rather than silently changing analyses.
message("Historical macro import completed.")
message("Samples: ", nrow(macro), "; sites: ", n_distinct(macro$site_id),
        "; taxa fields: ", length(taxa_cols))
message("Six-site study samples: ", nrow(macro_study),
        "; core historical samples: ", nrow(macro_core))
message("Processed data: ", processed_dir)
message("Inventories and QGIS export: ", output_dir)
