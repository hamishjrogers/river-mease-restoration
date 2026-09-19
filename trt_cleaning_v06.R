# Title: trt_cleaning_v06
# Author: HR
# Date: 2026-08-20
# Purpose: Harmonise phosphorus monitoring datasets originating from multiple
# Environment Agency monitoring programmes with the 2026 wetland restoration
# monitoring dataset into a single analytical dataset with consistent structure
# and reporting units.
#
# Inputs:
#   data/Mease_datasets_summary_v03.xlsx
#   data/wetland_2026_clean.rds
#   data/conversion_decision_2021.rds
#
# Outputs:
#   data/trt_clean_v06.csv
#   data/trt_clean_v06.rds
#   data/trt_site_id_dictionary_v06.csv
#   data/trt_site_id_dictionary_v06.rds

# ---- QC assumptions ----------------------------------------------------
# 1. Values reported below detection limit (BDL) are substituted using half
# the stated detection limit within each dataset independently.

# 2. Outlier detection and missingness diagnostics are conducted in separate
# quality-control workflows and are not retained as flags within the final
# analytical dataset.

# 3. Orthophosphate concentrations are standardised to mg/L as phosphorus (P)
# for comparability across monitoring programmes.

# 4. National Grid References (NGRs) are standardised and used to derive
# harmonised site identifiers across datasets.

# 5. Nitrite and ammonia were measured only during selected 2026 wetland
# restoration sampling rounds. These variables are therefore retained as NA
# for historic Environment Agency datasets.

# 6. The 2026 wetland restoration dataset contains only completed sampling
# observations. Placeholder rows are excluded upstream in
# wetland_restoration_cleaning_v03.R.
#
# 7. The reporting basis of 2021 orthophosphate is uncertain. The primary
# harmonisation adopts PO4 -> P conversion (/3.066), with the uncertainty
# retained explicitly in the diagnostic decision object and evaluated through
# sensitivity analysis.

# ---- 0: Libraries ------------------------------------------------------
# Core data handling and cleaning packages

library(readxl)
library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(readr)


# ---- 1: File path ------------------------------------------------------
# Path to master workbook containing all raw historic monitoring datasets

file_path <- "data/Mease_datasets_summary_v03.xlsx"


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


# Replace below-detection-limit values with DL/2.
# Character-based non-numeric values are converted to NA.
apply_bdl_half <- function(x, detection_limit) {
  
  x <- str_trim(as.character(x))
  
  case_when(
    str_detect(x, "<") ~ detection_limit / 2,
    x %in% c("", "NA", "N/A") ~ NA_real_,
    TRUE ~ suppressWarnings(as.numeric(x))
  )
}


# ---- 3: 2010 -----------------------------------------------------------
# Dataset contains orthophosphate already reported as P (mg/L)

trt_2010 <- read_excel(
  file_path,
  sheet = "2010 data"
) %>%
  
  # Standardise variable names
  clean_names() %>%
  
  # Rename key analytical variables
  rename(
    orthophosphate_as_p_mgl = `0180_orthophospht_mgl`,
    total_phosphorus_as_p_mgl = `0348_phosphorusp_mgl`
  ) %>%
  
  # Convert dates, handle BDL values, and standardise NGRs
  mutate(
    date_sampled = as.Date(date),
    
    orthophosphate_as_p_mgl =
      apply_bdl_half(orthophosphate_as_p_mgl, 0.02),
    
    total_phosphorus_as_p_mgl =
      apply_bdl_half(total_phosphorus_as_p_mgl, 0.02),
    
    year = 2010,
    ngr = standardise_ngr(ngr),
    
    source_programme = "EA_ROUTINE",
    site_type = "REFERENCE",
    
    site_id_original = NA_character_,
    site_name = NA_character_,
    survey_round = NA_integer_,
    time_sampled = NA,
    
    orthophosphate_as_po4_mgl = NA_real_,
    nitrite_no2_mgl = NA_real_,
    ammonia_nh3n_mgl = NA_real_,
    
    ph = NA_real_,
    dissolved_oxygen_mgl = NA_real_,
    temperature_c = NA_real_,
    oxygen_saturation_percent = NA_real_,
    pressure_hpa = NA_real_,
    
    easting = NA_real_,
    northing = NA_real_
  ) %>%
  
  # Retain harmonised analytical variables only
  select(
    site_id_original,
    site_name,
    tributary_name = ptname,
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
  )


# ---- 4: 2020 -----------------------------------------------------------
# Dataset contains orthophosphate already reported as P (mg/L)

# Raw sheet is stored in an inverted format requiring reshaping
trt_2020 <- read_excel(
  file_path,
  sheet = "WQ_2020",
  col_names = FALSE
)

# Assign variable column name
colnames(trt_2020)[1] <- "variable"

trt_2020 <- trt_2020 %>%
  
  # Reshape inverted structure into tidy tabular format
  pivot_longer(
    -variable,
    names_to = "sample_id",
    values_to = "value"
  ) %>%
  
  pivot_wider(
    names_from = variable,
    values_from = value
  ) %>%
  
  # Standardise variable names
  clean_names() %>%
  
  # Convert dates, analytical values, and standardise NGRs
  mutate(
    sample_taken_num = as.numeric(sample_taken),
    
    date_sampled =
      as.Date(sample_taken_num, origin = "1899-12-30"),
    
    orthophosphate_as_p_mgl =
      as.numeric(orthophosphate_reactive_as_p_mgl),
    
    total_phosphorus_as_p_mgl =
      as.numeric(phosphorus__total_as_p_mgl),
    
    year = 2020,
    ngr = standardise_ngr(grid_reference),
    
    source_programme = "EA_ROUTINE",
    site_type = "REFERENCE",
    
    site_id_original = NA_character_,
    site_name = NA_character_,
    survey_round = NA_integer_,
    time_sampled = NA,
    
    orthophosphate_as_po4_mgl = NA_real_,
    nitrite_no2_mgl = NA_real_,
    ammonia_nh3n_mgl = NA_real_,
    
    ph = NA_real_,
    dissolved_oxygen_mgl = NA_real_,
    temperature_c = NA_real_,
    oxygen_saturation_percent = NA_real_,
    pressure_hpa = NA_real_,
    
    easting = NA_real_,
    northing = NA_real_
  ) %>%
  
  # Retain harmonised analytical variables only
  select(
    site_id_original,
    site_name,
    tributary_name = water_body,
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
  )


# ---- 5: 2021 -----------------------------------------------------------
# The 2021 source sheet does not state whether orthophosphate is reported as
# PO4 or as P. The primary harmonisation workflow adopts PO4 -> P conversion
# using the decision object created by trt_2021_diagnostics_v03.R. This is an
# explicit analytical assumption and is evaluated separately by sensitivity
# analysis.

conversion_decision_2021 <-
  readRDS("data/conversion_decision_2021.rds")

# Confirm that the expected diagnostic decision object has been loaded
stopifnot(
  isTRUE(conversion_decision_2021$convert),
  identical(conversion_decision_2021$basis, "PO4_to_P"),
  identical(
    conversion_decision_2021$certainty,
    "uncertain_reporting_basis"
  ),
  isTRUE(all.equal(conversion_decision_2021$factor, 3.066))
)

trt_2021 <- read_excel(
  file_path,
  sheet = "Mease site data_2021"
) %>%
  
  # Standardise variable names
  clean_names() %>%
  
  # Convert dates, analytical values, and standardise NGRs
  mutate(
    orthophosphate_as_po4_mgl =
      suppressWarnings(as.numeric(orthophosphate_mean)),
    
    total_phosphorus_as_p_mgl =
      suppressWarnings(as.numeric(total_phosphate_mean)),
    
    date_sampled = as.Date(date),
    year = year(date_sampled),
    ngr = standardise_ngr(ngr),
    
    source_programme = "EA_ROUTINE",
    site_type = "REFERENCE",
    
    site_id_original = NA_character_,
    site_name = NA_character_,
    survey_round = NA_integer_,
    time_sampled = NA,
    
    nitrite_no2_mgl = NA_real_,
    ammonia_nh3n_mgl = NA_real_,
    
    ph = NA_real_,
    dissolved_oxygen_mgl = NA_real_,
    temperature_c = NA_real_,
    oxygen_saturation_percent = NA_real_,
    pressure_hpa = NA_real_,
    
    easting = NA_real_,
    northing = NA_real_,
    
    # Apply adopted PO4 -> P harmonisation assumption
    orthophosphate_as_p_mgl =
      orthophosphate_as_po4_mgl /
      conversion_decision_2021$factor
  ) %>%
  
  # Retain harmonised analytical variables only
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
  )


# ---- 6: 2025 -----------------------------------------------------------
# Dataset reports orthophosphate as PO4 requiring conversion to P

trt_2025 <- read_excel(
  file_path,
  sheet = "WQ_2025"
) %>%
  
  # Standardise variable names
  clean_names() %>%
  
  # Rename key analytical variable
  rename(
    orthophosphate_as_po4_mgl =
      orthophosphate_as_po_mgl
  ) %>%
  
  # Convert dates, handle BDL values, and standardise NGRs
  mutate(
    
    # Handle BDL values independently for each parameter
    orthophosphate_as_po4_mgl =
      apply_bdl_half(
        orthophosphate_as_po4_mgl,
        0.062
      ),
    
    total_phosphorus_as_p_mgl =
      apply_bdl_half(
        phosphorus_total_mgl,
        0.020
      ),
    
    # Convert PO4 to P
    orthophosphate_as_p_mgl =
      orthophosphate_as_po4_mgl / 3.066,
    
    date_sampled = as.Date(date_sampled),
    year = 2025,
    ngr = standardise_ngr(ngr),
    
    source_programme = "EA_ROUTINE",
    site_type = "REFERENCE",
    
    site_id_original = NA_character_,
    site_name = NA_character_,
    survey_round = NA_integer_,
    time_sampled = NA,
    
    nitrite_no2_mgl = NA_real_,
    ammonia_nh3n_mgl = NA_real_,
    
    ph = NA_real_,
    dissolved_oxygen_mgl = NA_real_,
    temperature_c = NA_real_,
    oxygen_saturation_percent = NA_real_,
    pressure_hpa = NA_real_,
    
    easting = NA_real_,
    northing = NA_real_
  ) %>%
  
  # Retain harmonised analytical variables only
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
  )


# ---- 7: 2026 wetland restoration monitoring ---------------------------
# Import separately cleaned restoration monitoring dataset created by
# wetland_restoration_cleaning_v03.R.

wetland_2026 <-
  readRDS("data/wetland_2026_clean.rds")


# ---- 7a: Validate 2026 structure ---------------------------------------

required_2026_columns <- c(
  "site_id_original",
  "site_name",
  "tributary_name",
  "ngr",
  "easting",
  "northing",
  "date_sampled",
  "time_sampled",
  "year",
  "source_programme",
  "site_type",
  "survey_round",
  "orthophosphate_as_po4_mgl",
  "orthophosphate_as_p_mgl",
  "total_phosphorus_as_p_mgl",
  "nitrite_no2_mgl",
  "ammonia_nh3n_mgl",
  "ph",
  "dissolved_oxygen_mgl",
  "temperature_c",
  "oxygen_saturation_percent",
  "pressure_hpa"
)

stopifnot(
  all(required_2026_columns %in% names(wetland_2026))
)

stopifnot(
  all(wetland_2026$source_programme ==
        "WETLAND_RESTORATION")
)

stopifnot(
  n_distinct(wetland_2026$site_id_original) == 6
)

stopifnot(
  !any(is.na(wetland_2026$date_sampled))
)

stopifnot(
  !any(is.na(wetland_2026$ngr))
)


# ---- 8: Validate harmonised column structure ---------------------------

harmonised_columns <- c(
  "site_id_original",
  "site_name",
  "tributary_name",
  "ngr",
  "easting",
  "northing",
  "date_sampled",
  "time_sampled",
  "year",
  "source_programme",
  "site_type",
  "survey_round",
  "orthophosphate_as_po4_mgl",
  "orthophosphate_as_p_mgl",
  "total_phosphorus_as_p_mgl",
  "nitrite_no2_mgl",
  "ammonia_nh3n_mgl",
  "ph",
  "dissolved_oxygen_mgl",
  "temperature_c",
  "oxygen_saturation_percent",
  "pressure_hpa"
)

stopifnot(
  identical(names(trt_2010), harmonised_columns),
  identical(names(trt_2020), harmonised_columns),
  identical(names(trt_2021), harmonised_columns),
  identical(names(trt_2025), harmonised_columns),
  identical(names(wetland_2026), harmonised_columns)
)


# ---- 9: Site ID dictionary --------------------------------------------
# Create harmonised site identifiers based on unique standardised NGRs.
# NGR remains the underlying stable definition of a physical monitoring site.

all_ngr <- bind_rows(
  trt_2010 %>% select(ngr),
  trt_2020 %>% select(ngr),
  trt_2021 %>% select(ngr),
  trt_2025 %>% select(ngr),
  wetland_2026 %>% select(ngr)
) %>%
  
  filter(
    !is.na(ngr),
    ngr != ""
  ) %>%
  
  distinct(ngr) %>%
  
  mutate(
    site_id = row_number()
  )


# ---- 9a: Site dictionary QA -------------------------------------------

stopifnot(
  !anyDuplicated(all_ngr$ngr)
)

stopifnot(
  !anyDuplicated(all_ngr$site_id)
)

stopifnot(
  !any(is.na(all_ngr$ngr))
)


# ---- 10: Join IDs and combine datasets --------------------------------
# Bind all harmonised datasets and append NGR-derived site identifiers.

trt_all_clean <- bind_rows(
  trt_2010,
  trt_2020,
  trt_2021,
  trt_2025,
  wetland_2026
) %>%
  
  left_join(
    all_ngr,
    by = "ngr"
  )


# Confirm that every retained record has received a harmonised site ID
stopifnot(
  !any(is.na(trt_all_clean$site_id))
)


# ---- 11: Tributary cleaning and classification -------------------------

trt_all_clean <- trt_all_clean %>%
  
  mutate(
    
    # Standardise tributary names first
    tributary_name_clean = tributary_name %>%
      str_to_upper() %>%
      str_replace_all("\\b\\d+\\b", "") %>%
      str_squish(),
    
    # Assign tributary family
    tributary_family = case_when(
      
      # Mease system: highest priority
      str_detect(
        tributary_name_clean,
        "CROXALL"
      ) ~ "MEASE - CROXALL",
      
      str_detect(
        tributary_name_clean,
        "CLIFTON CAMPVILLE"
      ) ~ "MEASE - CLIFTON CAMPVILLE",
      
      str_detect(
        tributary_name_clean,
        "STRETTON BRIDGE"
      ) ~ "MEASE - STRETTON BRIDGE",
      
      str_detect(
        tributary_name_clean,
        "BIRDSHILL"
      ) ~ "MEASE - BIRDSHILL ROAD",
      
      str_detect(
        tributary_name_clean,
        "SNARESTONE"
      ) ~ "MEASE - SNARESTONE STW",
      
      str_detect(
        tributary_name_clean,
        "UPPER MEASE"
      ) ~ "MEASE - UPPER MEASE",
      
      str_detect(
        tributary_name_clean,
        "MEASE-GILWISKAW CONFLUENCE"
      ) ~ "MEASE - CONFLUENCE",
      
      str_detect(
        tributary_name_clean,
        "^RIVER MEASE$"
      ) ~ "MEASE - MAIN STEM",
      
      str_detect(
        tributary_name_clean,
        "^MEASE$"
      ) ~ "MEASE - MAIN STEM",
      
      str_detect(
        tributary_name_clean,
        "MEASE"
      ) ~ "MEASE - OTHER",
      
      # Tributaries
      str_detect(
        tributary_name_clean,
        "APPLEBY"
      ) ~ "APPLEBY BROOK",
      
      str_detect(
        tributary_name_clean,
        "CHILCOTE"
      ) ~ "CHILCOTE BROOK",
      
      str_detect(
        tributary_name_clean,
        "GILWISKAW"
      ) ~ "GILWISKAW BROOK",
      
      str_detect(
        tributary_name_clean,
        "HOOBOROUGH"
      ) ~ "HOOBOROUGH BROOK",
      
      str_detect(
        tributary_name_clean,
        "HARLASTON"
      ) ~ "HARLASTON BROOK",
      
      str_detect(
        tributary_name_clean,
        "KES"
      ) ~ "KES BROOK",
      
      str_detect(
        tributary_name_clean,
        "PESSAL"
      ) ~ "PESSAL BROOK",
      
      str_detect(
        tributary_name_clean,
        "SALTERSFORD"
      ) ~ "SALTERSFORD BROOK",
      
      str_detect(
        tributary_name_clean,
        "SHELL"
      ) ~ "SHELL BROOK",
      
      str_detect(
        tributary_name_clean,
        "SEAL"
      ) ~ "SEAL BROOK",
      
      str_detect(
        tributary_name_clean,
        "WEST"
      ) ~ "WEST BROOK",
      
      str_detect(
        tributary_name_clean,
        "HARLASTAN"
      ) ~ "HARLASTON BROOK",
      
      str_detect(
        tributary_name_clean,
        "ALTON GRANGE"
      ) ~ "GILWISKAW BROOK",
      
      str_detect(
        tributary_name_clean,
        "FLAGSTAFF"
      ) ~ "GILWISKAW BROOK",
      
      str_detect(
        tributary_name_clean,
        "NORMANTON LE HEATH"
      ) ~ "GILWISKAW BROOK",
      
      str_detect(
        tributary_name_clean,
        "PACKINGTON"
      ) ~ "GILWISKAW BROOK",
      
      str_detect(
        tributary_name_clean,
        "SMISBY"
      ) ~ "GILWISKAW BROOK",
      
      str_detect(
        tributary_name_clean,
        "ASHBY CANAL"
      ) ~ "ASHBY CANAL",
      
      TRUE ~ NA_character_
    )
  )


# ---- 12: Final QA ------------------------------------------------------

summary(trt_all_clean)

nrow(trt_all_clean)

n_distinct(trt_all_clean$site_id)

table(
  trt_all_clean$year,
  useNA = "ifany"
)

table(
  trt_all_clean$source_programme,
  useNA = "ifany"
)

table(
  trt_all_clean$survey_round,
  useNA = "ifany"
)


# Confirm expected source programmes
stopifnot(
  setequal(
    unique(trt_all_clean$source_programme),
    c(
      "EA_ROUTINE",
      "WETLAND_RESTORATION"
    )
  )
)


# Confirm that only the six intended contemporary sites are represented
restoration_site_check <- trt_all_clean %>%
  filter(
    source_programme ==
      "WETLAND_RESTORATION"
  ) %>%
  distinct(
    site_id_original,
    site_name,
    ngr,
    site_id
  ) %>%
  arrange(site_id_original)

restoration_site_check

stopifnot(
  nrow(restoration_site_check) == 6
)


# Confirm that every completed contemporary round contains six sites
restoration_round_check <- trt_all_clean %>%
  filter(
    source_programme ==
      "WETLAND_RESTORATION"
  ) %>%
  count(
    survey_round,
    date_sampled,
    name = "n_sites"
  )

restoration_round_check

stopifnot(
  all(restoration_round_check$n_sites == 6)
)


# Review missingness by source programme
missingness_by_source <- trt_all_clean %>%
  group_by(source_programme) %>%
  summarise(
    across(
      everything(),
      ~ sum(is.na(.x))
    ),
    .groups = "drop"
  )

missingness_by_source


# Review tributary names not assigned to a tributary family
unclassified_tributaries <- trt_all_clean %>%
  filter(
    is.na(tributary_family)
  ) %>%
  distinct(
    tributary_name,
    tributary_name_clean
  ) %>%
  arrange(tributary_name_clean)

unclassified_tributaries


# ---- 13: Save ----------------------------------------------------------

write.csv(
  trt_all_clean,
  "data/trt_clean_v06.csv",
  row.names = FALSE
)

saveRDS(
  trt_all_clean,
  "data/trt_clean_v06.rds"
)


# Save NGR-derived site dictionary for transparent downstream reference
write.csv(
  all_ngr,
  "data/trt_site_id_dictionary_v06.csv",
  row.names = FALSE
)

saveRDS(
  all_ngr,
  "data/trt_site_id_dictionary_v06.rds"
)

# ---- 14: Historical orthophosphate censoring QC ------------------------

# Quantify raw below-detection-limit (BDL) orthophosphate observations used in
# the historical analysis. Explicit "<" values occur in the 2010 and 2025
# source datasets; 2020 and 2021 are checked separately for censored or
# non-numeric values.

# 2010 raw orthophosphate
bdl_2010 <- read_excel(
  file_path,
  sheet = "2010 data"
) %>%
  clean_names() %>%
  summarise(
    n_ortho = sum(!is.na(`0180_orthophospht_mgl`)),
    n_bdl = sum(
      str_detect(
        str_trim(as.character(`0180_orthophospht_mgl`)),
        "<"
      ),
      na.rm = TRUE
    )
  )

# 2025 raw orthophosphate
bdl_2025 <- read_excel(
  file_path,
  sheet = "WQ_2025"
) %>%
  clean_names() %>%
  summarise(
    n_ortho = sum(!is.na(orthophosphate_as_po_mgl)),
    n_bdl = sum(
      str_detect(
        str_trim(as.character(orthophosphate_as_po_mgl)),
        "<"
      ),
      na.rm = TRUE
    )
  )

bdl_summary <- bind_rows(
  `2010` = bdl_2010,
  `2025` = bdl_2025,
  .id = "dataset"
) %>%
  mutate(
    percent_bdl = 100 * n_bdl / n_ortho
  )

# Reconstruct raw 2020 structure and check for censored / non-numeric values
raw_2020 <- read_excel(
  file_path,
  sheet = "WQ_2020",
  col_names = FALSE
)

colnames(raw_2020)[1] <- "variable"

raw_2020_tidy <- raw_2020 %>%
  pivot_longer(
    -variable,
    names_to = "sample_id",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = variable,
    values_from = value
  ) %>%
  clean_names()

check_2020 <- raw_2020_tidy %>%
  summarise(
    n_nonmissing = sum(!is.na(orthophosphate_reactive_as_p_mgl)),
    n_less_than = sum(
      str_detect(
        str_trim(as.character(orthophosphate_reactive_as_p_mgl)),
        "<"
      ),
      na.rm = TRUE
    ),
    n_nonnumeric = sum(
      !is.na(orthophosphate_reactive_as_p_mgl) &
        is.na(
          suppressWarnings(
            as.numeric(orthophosphate_reactive_as_p_mgl)
          )
        )
    )
  )

# 2021: verify that non-numeric entries are missing-value codes rather than BDL
raw_2021 <- read_excel(
  file_path,
  sheet = "Mease site data_2021"
) %>%
  clean_names()

check_2021 <- raw_2021 %>%
  summarise(
    n_nonmissing = sum(!is.na(orthophosphate_mean)),
    n_less_than = sum(
      str_detect(
        str_trim(as.character(orthophosphate_mean)),
        "<"
      ),
      na.rm = TRUE
    ),
    n_nonnumeric = sum(
      !is.na(orthophosphate_mean) &
        is.na(
          suppressWarnings(
            as.numeric(orthophosphate_mean)
          )
        )
    )
  )

non_numeric_2021 <- raw_2021 %>%
  filter(
    !is.na(orthophosphate_mean),
    is.na(
      suppressWarnings(
        as.numeric(orthophosphate_mean)
      )
    )
  ) %>%
  count(
    orthophosphate_mean,
    sort = TRUE
  )

# Confirm expected raw-data QC findings
stopifnot(
  bdl_2010$n_bdl == 26,
  bdl_2025$n_bdl == 21,
  check_2020$n_less_than == 0,
  check_2020$n_nonnumeric == 0,
  check_2021$n_less_than == 0,
  check_2021$n_nonnumeric == 22,
  nrow(non_numeric_2021) == 1,
  non_numeric_2021$orthophosphate_mean[1] == "n/a",
  non_numeric_2021$n[1] == 22
)

# Historical analytical dataset excludes contemporary restoration monitoring.
historic_ortho_n <- trt_all_clean %>%
  filter(
    source_programme == "EA_ROUTINE",
    !is.na(orthophosphate_as_p_mgl)
  ) %>%
  nrow()

bdl_historic_analysis <- tibble(
  n_bdl = sum(bdl_summary$n_bdl),
  n_historic_orthophosphate = historic_ortho_n,
  percent_of_historic_dataset =
    100 * n_bdl / n_historic_orthophosphate
)

# Confirm dissertation-reported historical sample count and BDL proportion
stopifnot(
  historic_ortho_n == 633,
  bdl_historic_analysis$n_bdl == 47
)


# ---- 15: Final QC summary ----------------------------------------------

cat("\nTRT harmonised cleaning complete\n")
cat("--------------------------------\n")
cat("Total rows:", nrow(trt_all_clean), "\n")
cat(
  "Historic usable orthophosphate observations:",
  historic_ortho_n,
  "\n"
)
cat(
  "Contemporary restoration observations:",
  sum(
    trt_all_clean$source_programme ==
      "WETLAND_RESTORATION"
  ),
  "\n"
)
cat(
  "Unique harmonised sites:",
  n_distinct(trt_all_clean$site_id),
  "\n"
)
cat(
  "Historic BDL orthophosphate observations:",
  bdl_historic_analysis$n_bdl,
  sprintf(
    "(%.2f%%)\n",
    bdl_historic_analysis$percent_of_historic_dataset
  )
)
cat(
  "2021 reporting-basis certainty:",
  conversion_decision_2021$certainty,
  "\n\n"
)

print(
  table(
    trt_all_clean$year,
    useNA = "ifany"
  )
)

print(
  table(
    trt_all_clean$source_programme,
    useNA = "ifany"
  )
)

print(restoration_site_check)
print(restoration_round_check)
print(bdl_summary)
print(check_2020)
print(check_2021)
print(non_numeric_2021)
print(bdl_historic_analysis)
