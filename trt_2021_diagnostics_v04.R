# ---- Header ------------------------------------------------------------
# Title: trt_2021_diagnostics_v04
# Author: HR
# Date: 2026-08-20
# Purpose: Evaluate the reporting-basis uncertainty in the 2021 orthophosphate
# dataset before harmonisation with other River Mease phosphorus datasets.
#
# Input:
#   data/Mease_datasets_summary_v03.xlsx
#   sheet: "Mease site data_2021"
#
# Output:
#   data/conversion_decision_2021.rds
#
# The saved decision object is read by the downstream historic phosphorus
# cleaning workflow.

# ---- Diagnostic assumptions -------------------------------------------
# 1. The 2021 source sheet labels the relevant variables as
#    "Orthophosphate (Mean)" and "Total Phosphate (Mean)" but does not state
#    whether orthophosphate is reported as PO4 or as phosphorus (P).
#
# 2. Total phosphorus values are treated as phosphorus (P; mg/L) in the
#    downstream harmonisation workflow.
#
# 3. Two candidate interpretations of the raw 2021 orthophosphate values are
#    considered:
#       a) values already reported as P
#       b) values reported as PO4 and converted to P using PO4 / 3.066
#
# 4. Ratio analysis, visual inspection and linear regression are retained as
#    supporting diagnostics only. These diagnostics cannot determine the
#    reporting basis conclusively because scaling orthophosphate by a constant
#    does not alter the strength of its relationship with total phosphorus.
#
# 5. For the primary harmonisation workflow, the raw 2021 orthophosphate values
#    are treated as PO4 and converted to P using / 3.066. This is an explicit
#    analytical assumption rather than a proven unit identification.
#
# 6. The effect of this assumption is evaluated separately through sensitivity
#    analysis using alternative unit interpretations.
#
# 7. Non-numeric source entries (e.g. "n/a") are converted to NA and excluded
#    automatically from numerical diagnostics.

# ---- 0: Libraries ------------------------------------------------------

library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)

# ---- 1: File path ------------------------------------------------------

file_path <- "data/Mease_datasets_summary_v03.xlsx"

# ---- 2: Helper function ------------------------------------------------

clean_names <- function(df) {
  df %>%
    rename_with(
      ~ .x %>%
        str_to_lower() %>%
        str_replace_all("\\s+", "_") %>%
        str_replace_all("[^a-z0-9_]", "")
    )
}

# ---- 3: Import and prepare 2021 dataset --------------------------------

trt_2021 <- read_excel(
  file_path,
  sheet = "Mease site data_2021"
) %>%
  clean_names() %>%
  mutate(
    orthophosphate_raw_mgl =
      suppressWarnings(as.numeric(orthophosphate_mean)),
    total_phosphorus_as_p_mgl =
      suppressWarnings(as.numeric(total_phosphate_mean))
  )

stopifnot(
  all(
    c(
      "orthophosphate_raw_mgl",
      "total_phosphorus_as_p_mgl"
    ) %in% names(trt_2021)
  )
)

n_paired <- trt_2021 %>%
  filter(
    !is.na(orthophosphate_raw_mgl),
    !is.na(total_phosphorus_as_p_mgl)
  ) %>%
  nrow()

stopifnot(n_paired > 0)

# ---- 4: Candidate reporting scenarios ----------------------------------

trt_2021_diag <- trt_2021 %>%
  mutate(
    ortho_as_p_no_conversion =
      orthophosphate_raw_mgl,

    ortho_as_p_if_po4 =
      orthophosphate_raw_mgl / 3.066,

    ratio_no_conversion =
      if_else(
        total_phosphorus_as_p_mgl > 0,
        ortho_as_p_no_conversion / total_phosphorus_as_p_mgl,
        NA_real_
      ),

    ratio_if_po4 =
      if_else(
        total_phosphorus_as_p_mgl > 0,
        ortho_as_p_if_po4 / total_phosphorus_as_p_mgl,
        NA_real_
      )
  )

# ---- 5: Ratio diagnostics ----------------------------------------------

ratio_summary <- trt_2021_diag %>%
  summarise(
    n_paired = sum(
      !is.na(orthophosphate_raw_mgl) &
        !is.na(total_phosphorus_as_p_mgl)
    ),

    median_ratio_no_conversion =
      median(ratio_no_conversion, na.rm = TRUE),

    median_ratio_if_po4 =
      median(ratio_if_po4, na.rm = TRUE),

    proportion_ortho_gt_tp_no_conversion =
      mean(
        ortho_as_p_no_conversion > total_phosphorus_as_p_mgl,
        na.rm = TRUE
      ),

    proportion_ortho_gt_tp_if_po4 =
      mean(
        ortho_as_p_if_po4 > total_phosphorus_as_p_mgl,
        na.rm = TRUE
      )
  )

ratio_summary

# ---- 6: Visual inspection ----------------------------------------------

diagnostic_plot <- trt_2021_diag %>%
  select(
    total_phosphorus_as_p_mgl,
    ortho_as_p_no_conversion,
    ortho_as_p_if_po4
  ) %>%
  pivot_longer(
    cols = c(
      ortho_as_p_no_conversion,
      ortho_as_p_if_po4
    ),
    names_to = "scenario",
    values_to = "orthophosphate_as_p_mgl"
  ) %>%
  mutate(
    scenario = recode(
      scenario,
      ortho_as_p_no_conversion =
        "No conversion (raw treated as P)",
      ortho_as_p_if_po4 =
        "PO4 converted to P (/3.066)"
    )
  ) %>%
  ggplot(
    aes(
      x = total_phosphorus_as_p_mgl,
      y = orthophosphate_as_p_mgl
    )
  ) +
  geom_point(alpha = 0.5) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE
  ) +
  facet_wrap(~ scenario) +
  labs(
    x = "Total phosphorus (mg/L as P)",
    y = "Orthophosphate (mg/L as P)"
  )

diagnostic_plot

# ---- 7: Regression diagnostics -----------------------------------------

lm_no_conversion <- lm(
  ortho_as_p_no_conversion ~ total_phosphorus_as_p_mgl,
  data = trt_2021_diag
)

lm_if_po4 <- lm(
  ortho_as_p_if_po4 ~ total_phosphorus_as_p_mgl,
  data = trt_2021_diag
)

regression_summary <- tibble(
  scenario = c(
    "No conversion (raw treated as P)",
    "PO4 converted to P (/3.066)"
  ),
  intercept = c(
    unname(coef(lm_no_conversion)[1]),
    unname(coef(lm_if_po4)[1])
  ),
  slope = c(
    unname(coef(lm_no_conversion)[2]),
    unname(coef(lm_if_po4)[2])
  ),
  r_squared = c(
    summary(lm_no_conversion)$r.squared,
    summary(lm_if_po4)$r.squared
  )
)

regression_summary

# ---- 8: Harmonisation decision -----------------------------------------

# The available metadata and diagnostic relationships do not determine the
# original reporting basis conclusively.
#
# For the primary harmonisation workflow, the raw 2021 orthophosphate values
# are treated as PO4 and converted to phosphorus equivalents using / 3.066.
# This is an explicit analytical assumption. Alternative interpretations are
# evaluated separately through sensitivity analysis.

conversion_decision_2021 <- list(
  convert = TRUE,
  factor = 3.066,
  basis = "PO4_to_P",
  certainty = "uncertain_reporting_basis",

  evidence = list(
    n_paired = n_paired,
    ratio_summary = ratio_summary,
    regression_summary = regression_summary,
    r2_raw_relationship =
      summary(lm_no_conversion)$r.squared
  )
)

stopifnot(
  isTRUE(conversion_decision_2021$convert),
  identical(conversion_decision_2021$basis, "PO4_to_P"),
  identical(
    conversion_decision_2021$certainty,
    "uncertain_reporting_basis"
  ),
  isTRUE(all.equal(conversion_decision_2021$factor, 3.066))
)

# ---- 9: Save ------------------------------------------------------------

saveRDS(
  conversion_decision_2021,
  "data/conversion_decision_2021.rds"
)

# ---- 10: Final QC summary ----------------------------------------------

cat("\n2021 orthophosphate diagnostic complete\n")
cat("---------------------------------------\n")
cat("Paired orthophosphate/TP observations:", n_paired, "\n")
cat(
  "Adopted primary reporting basis:",
  conversion_decision_2021$basis,
  "\n"
)
cat(
  "Decision certainty:",
  conversion_decision_2021$certainty,
  "\n"
)
cat("Conversion factor:", conversion_decision_2021$factor, "\n\n")

print(ratio_summary)
print(regression_summary)
