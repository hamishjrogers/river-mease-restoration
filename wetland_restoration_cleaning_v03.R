# ---- Header ------------------------------------------------------------
# Title: wetland_restoration_cleaning_v03
# Author: HR
# Date: 2026-08-20
# Purpose: Clean and harmonise the final 2026 wetland-restoration water-quality
# dataset for integration with historic Environment Agency phosphorus monitoring
# data. The final analytical dataset contains six sites sampled across three rounds.
#
# Input:
#   data/2026_measewqsampling.csv
#
# Outputs:
#   data/wetland_2026_clean.csv
#   data/wetland_2026_clean.rds
#   data/wetland_2026_placeholders.csv

# ---- QC assumptions ----------------------------------------------------
# 1. Orthophosphate concentrations reported as PO4 (mg/L) are accompanied
# by concentrations converted to phosphorus as P (mg/L).

# 2. Field measurements are spot measurements collected contemporaneously
# with nutrient samples.

# 3. National Grid References (NGRs) are standardised to ensure consistency
# across monitoring programmes.

# 4. Missing or non-numeric analytical values are converted to NA.

# 5. No below-detection-limit (BDL) substitution is applied to the 2026 dataset.
# Raw PO4 and P fields are checked for explicit "<" censoring before numeric
# conversion so that censored observations cannot be silently converted to NA.

# 6. Rows without a sampling date are treated as placeholders, retained separately
# for transparent QA, and excluded from the analytical dataset.

# 7. Sampling-round numbers are derived from the chronological sequence
# of completed sampling dates rather than being manually assigned.

# ---- 0: Libraries ------------------------------------------------------

library(readr)
library(dplyr)
library(stringr)
library(lubridate)


# ---- 1: File path ------------------------------------------------------

file_path <- "data/2026_measewqsampling.csv"


# ---- 2: Helper functions -----------------------------------------------

# Standardise column names
clean_names <- function(df) {
  df %>%
    rename_with(
      ~ .x %>%
        str_to_lower() %>%
        str_replace_all("\\s+", "_") %>%
        str_replace_all("[^a-z0-9_]", "")
    )
}


# Standardise National Grid References
standardise_ngr <- function(x) {
  x %>%
    str_to_upper() %>%
    str_replace_all("\\s+", "") %>%
    str_trim()
}


# Convert character or numeric fields safely to numeric
as_numeric_safe <- function(x) {
  suppressWarnings(as.numeric(x))
}


# ---- 3: Import raw dataset ---------------------------------------------

wetland_2026_raw <- read_csv(
  file_path,
  show_col_types = FALSE,
  na = c("", "NA", "N/A")
) %>%
  
  # Standardise raw column names
  clean_names()


# Inspect imported structure
glimpse(wetland_2026_raw)

names(wetland_2026_raw)


# ---- 3a: Validate raw structure ----------------------------------------

required_raw_columns <- c(
  "site_id",
  "descriptive_site_id",
  "date",
  "ngr",
  "easting",
  "northing",
  "po4_mgl",
  "p_mgl",
  "ph",
  "do_concentration_mgl",
  "temperature_c",
  "oxygen_saturation_",
  "pressure_hpa",
  "nitrite_n02_mgl",
  "ammonia_nh3n_mgl"
)

stopifnot(
  all(required_raw_columns %in% names(wetland_2026_raw))
)

# Confirm that only the intended six monitoring sites are present
expected_sites <- c(
  "Upstream Gilwiskaw",
  "Wetland Gilwiskaw",
  "Upstream Mease",
  "Wetland Mease",
  "Confluence",
  "Downstream Mease"
)

stopifnot(
  setequal(
    unique(wetland_2026_raw$descriptive_site_id),
    expected_sites
  )
)

# Confirm that all site codes are present
expected_site_codes <- c(
  "UG",
  "WG",
  "UM",
  "WM",
  "C",
  "DM"
)

stopifnot(
  setequal(
    unique(wetland_2026_raw$site_id),
    expected_site_codes
  )
)


# Confirm that 2026 orthophosphate fields do not contain explicit BDL censoring.
# If future raw data include values such as "<0.05", stop here rather than
# allowing numeric conversion to turn them silently into NA.
bdl_2026_check <- wetland_2026_raw %>%
  summarise(
    po4_bdl = sum(str_detect(str_trim(as.character(po4_mgl)), "<"), na.rm = TRUE),
    p_bdl = sum(str_detect(str_trim(as.character(p_mgl)), "<"), na.rm = TRUE)
  )

print(bdl_2026_check)

stopifnot(
  bdl_2026_check$po4_bdl == 0,
  bdl_2026_check$p_bdl == 0
)


# ---- 4: Identify placeholder rows --------------------------------------

# Placeholder rows contain a site identifier but do not yet contain a
# sampling date. They are retained separately for QA and excluded from the
# analytical dataset until field data have been entered.

wetland_2026_placeholders <- wetland_2026_raw %>%
  filter(is.na(date))

wetland_2026_completed_raw <- wetland_2026_raw %>%
  filter(!is.na(date))


# QA: inspect completed and placeholder row counts
nrow(wetland_2026_completed_raw)

nrow(wetland_2026_placeholders)

wetland_2026_placeholders %>%
  count(descriptive_site_id)


# Every represented sampling date should contain exactly six sites
completed_date_counts <- wetland_2026_completed_raw %>%
  count(date, name = "n_sites")

completed_date_counts

stopifnot(
  all(completed_date_counts$n_sites == 6)
)


# ---- 5: Clean completed observations -----------------------------------

wetland_2026 <- wetland_2026_completed_raw %>%
  
  # Rename variables to harmonised structure
  rename(
    site_id_original = site_id,
    site_name = descriptive_site_id,
    orthophosphate_as_po4_mgl = po4_mgl,
    orthophosphate_as_p_mgl = p_mgl,
    dissolved_oxygen_mgl = do_concentration_mgl,
    oxygen_saturation_percent = oxygen_saturation_,
    nitrite_no2_mgl = nitrite_n02_mgl,
    ammonia_nh3n_mgl = ammonia_nh3n_mgl
  ) %>%
  
  # Convert variables and standardise formats
  mutate(
    
    # Date conversion
    date_sampled = as.Date(date),
    
    # No sampling-time variable is currently present in the raw CSV
    time_sampled = NA,
    
    # Year extraction
    year = year(date_sampled),
    
    # Standardise site identifiers and names
    site_id_original = str_to_upper(str_trim(site_id_original)),
    site_name = str_squish(site_name),
    
    # Standardise National Grid Reference
    ngr = standardise_ngr(ngr),
    
    # Convert coordinates to numeric
    easting = as_numeric_safe(easting),
    northing = as_numeric_safe(northing),
    
    # Convert nutrient variables to numeric
    orthophosphate_as_po4_mgl =
      as_numeric_safe(orthophosphate_as_po4_mgl),
    
    orthophosphate_as_p_mgl =
      as_numeric_safe(orthophosphate_as_p_mgl),
    
    nitrite_no2_mgl =
      as_numeric_safe(nitrite_no2_mgl),
    
    ammonia_nh3n_mgl =
      as_numeric_safe(ammonia_nh3n_mgl),
    
    # Convert field variables to numeric
    ph = as_numeric_safe(ph),
    
    dissolved_oxygen_mgl =
      as_numeric_safe(dissolved_oxygen_mgl),
    
    temperature_c =
      as_numeric_safe(temperature_c),
    
    oxygen_saturation_percent =
      as_numeric_safe(oxygen_saturation_percent),
    
    pressure_hpa =
      as_numeric_safe(pressure_hpa),
    
    # Derive sampling round from chronological sampling date
    survey_round =
      dense_rank(date_sampled),
    
    # Create broad watercourse classification
    tributary_name = case_when(
      site_name %in% c(
        "Upstream Gilwiskaw",
        "Wetland Gilwiskaw"
      ) ~ "Gilwiskaw Brook",
      
      site_name %in% c(
        "Upstream Mease",
        "Wetland Mease",
        "Downstream Mease"
      ) ~ "River Mease",
      
      site_name == "Confluence" ~
        "Mease-Gilwiskaw Confluence",
      
      TRUE ~ NA_character_
    ),
    
    # Metadata fields
    source_programme = "WETLAND_RESTORATION",
    site_type = "RESTORATION_MONITORING",
    
    # Historic dataset contains total phosphorus, but it was not measured
    # during the contemporary restoration survey
    total_phosphorus_as_p_mgl = NA_real_
  ) %>%
  
  # Retain harmonised analytical variables
  select(
    site_id_original,
    site_name,
    tributary_name,
    ngr,
    easting,
    northing,
    date_sampled,
    time_sampled,
    year,
    source_programme,
    site_type,
    survey_round,
    orthophosphate_as_po4_mgl,
    orthophosphate_as_p_mgl,
    total_phosphorus_as_p_mgl,
    nitrite_no2_mgl,
    ammonia_nh3n_mgl,
    ph,
    dissolved_oxygen_mgl,
    temperature_c,
    oxygen_saturation_percent,
    pressure_hpa
  ) %>%
  
  # Apply intended longitudinal site order
  mutate(
    site_name = factor(
      site_name,
      levels = c(
        "Upstream Gilwiskaw",
        "Wetland Gilwiskaw",
        "Upstream Mease",
        "Wetland Mease",
        "Confluence",
        "Downstream Mease"
      )
    )
  ) %>%
  
  arrange(
    survey_round,
    site_name
  )


# ---- 6: Cleaned-data validation ----------------------------------------

glimpse(wetland_2026)

summary(wetland_2026)


# Confirm exactly six sites are represented
stopifnot(
  n_distinct(wetland_2026$site_id_original) == 6
)

stopifnot(
  n_distinct(wetland_2026$site_name) == 6
)


# Confirm the final 2026 design: 18 observations across three sampling rounds
stopifnot(
  nrow(wetland_2026) == 18,
  n_distinct(wetland_2026$survey_round) == 3
)


# Confirm that every completed round contains six sites
round_counts <- wetland_2026 %>%
  count(
    survey_round,
    date_sampled,
    name = "n_sites"
  )

round_counts

stopifnot(
  all(round_counts$n_sites == 6)
)


# Confirm that each site occurs only once per sampling date
duplicate_site_dates <- wetland_2026 %>%
  count(
    site_id_original,
    date_sampled,
    name = "n"
  ) %>%
  filter(n > 1)

duplicate_site_dates

stopifnot(
  nrow(duplicate_site_dates) == 0
)


# Confirm that site codes consistently match site names
site_dictionary_check <- wetland_2026 %>%
  distinct(
    site_id_original,
    site_name
  )

site_dictionary_check

stopifnot(
  nrow(site_dictionary_check) == 6
)


# Confirm that spatial identifiers are complete before checking within-site consistency
stopifnot(
  !any(is.na(wetland_2026$ngr)),
  !any(is.na(wetland_2026$easting)),
  !any(is.na(wetland_2026$northing))
)


# Confirm that each site has one consistent NGR and coordinate pair
site_ngr_check <- wetland_2026 %>%
  group_by(
    site_id_original,
    site_name
  ) %>%
  summarise(
    n_ngr = n_distinct(ngr),
    n_easting = n_distinct(easting),
    n_northing = n_distinct(northing),
    .groups = "drop"
  )

site_ngr_check

stopifnot(
  all(site_ngr_check$n_ngr == 1),
  all(site_ngr_check$n_easting == 1),
  all(site_ngr_check$n_northing == 1)
)


# Confirm that PO4-to-P values are approximately consistent with the
# molecular conversion factor used elsewhere in the workflow
po4_conversion_check <- wetland_2026 %>%
  mutate(
    p_calculated = orthophosphate_as_po4_mgl / 3.066,
    conversion_difference =
      abs(orthophosphate_as_p_mgl - p_calculated)
  ) %>%
  select(
    site_name,
    date_sampled,
    orthophosphate_as_po4_mgl,
    orthophosphate_as_p_mgl,
    p_calculated,
    conversion_difference
  )

po4_conversion_check

stopifnot(
  all(
    po4_conversion_check$conversion_difference <= 0.002,
    na.rm = TRUE
  )
)


# Review missingness
colSums(is.na(wetland_2026))


# Review completed sampling rounds
table(wetland_2026$survey_round)

table(wetland_2026$date_sampled)


# ---- 7: Save -----------------------------------------------------------

write.csv(
  wetland_2026,
  "data/wetland_2026_clean.csv",
  row.names = FALSE
)

saveRDS(
  wetland_2026,
  "data/wetland_2026_clean.rds"
)


# Save placeholder rows separately for transparent QA
write.csv(
  wetland_2026_placeholders,
  "data/wetland_2026_placeholders.csv",
  row.names = FALSE
)

# ---- 8: Final reproducibility summary ----------------------------------

# Print compact end-of-script QA outputs without relying on hard-coded console
# comments from earlier versions of the dataset.
final_qc <- list(
  n_observations = nrow(wetland_2026),
  n_sites = n_distinct(wetland_2026$site_id_original),
  n_rounds = n_distinct(wetland_2026$survey_round),
  n_placeholders = nrow(wetland_2026_placeholders),
  round_counts = round_counts,
  bdl_check = bdl_2026_check,
  site_dictionary = site_dictionary_check,
  site_ngr_check = site_ngr_check,
  missing_values = colSums(is.na(wetland_2026))
)

final_qc

wfd-2000-6055-ec
