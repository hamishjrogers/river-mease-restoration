# ============================================================================
# River Mease | Historical macroinvertebrates
# Script 04: Taxonomic audit and provisional EPT richness
# File: macro_analysis_step03_v02.R
# ============================================================================
# PURPOSE
#   Audit the cleaned taxon names, identify provisional Ephemeroptera,
#   Plecoptera and Trichoptera (EPT) columns, calculate sample-level EPT
#   richness and describe variation across six study and three long-term sites.
#
# RUN ORDER / WORKING DIRECTORY
#   Run after macro_historical_v03.R, macro_analysis_step01_v02.R and
#   macro_analysis_step02_v02.R, from the EA_data project root.
#
# INPUTS
#   data/processed/macro_metrics_v01.rds
#   data/processed/macro_study_metrics_v01.rds
#   data/processed/macro_longterm_v01.rds
#   data/processed/macro_taxa_cols_v01.rds
#
# OUTPUTS
#   Retains original step03 CSV/PNG and *_ept_v01.rds filenames.
#   Adds CSV summaries of exploratory tests and an EPT candidate audit.
#
# LIMITATIONS
#   EPT names are selected by provisional keyword matching, NOT an externally
#   verified taxonomic lookup. Inspect step03_provisional_ept_taxa.csv and
#   step03_ept_candidate_audit.csv before interpreting EPT richness.
#   Excludes the known false-positive epoicocladius_ephemerae.
#   Counts source taxon columns present (>0); mixed taxonomic resolutions may
#   count overlapping higher/lower taxa. Source NA is treated as not present
#   for this metric, as in the original script.
#   Site tests are exploratory: uneven effort and repeated sampling violate
#   simple independent-sample assumptions. This script does NOT calculate BMWP
#   or ASPT: those require a verified family mapping and score lookup.
#   Revised from macro_analysis_step03_v01.R; validate against the source
#   analysis before treating the outputs as final.
# ============================================================================

# ---- 0. Packages and input checks ---------------------------------------
library(tidyverse)

processed_dir <- "data/processed"
output_dir <- "outputs/macroinvertebrates"
out <- function(name) file.path(output_dir, name)
inputs <- c("macro_metrics_v01.rds", "macro_study_metrics_v01.rds",
            "macro_longterm_v01.rds", "macro_taxa_cols_v01.rds")
missing_inputs <- inputs[!file.exists(file.path(processed_dir, inputs))]
if (length(missing_inputs)) stop(
  "Missing processed inputs: ", paste(missing_inputs, collapse = ", "),
  "\nRun Scripts 1–3 from the EA_data project root."
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

macro_metrics <- readRDS(file.path(processed_dir, inputs[1]))
macro_study_metrics <- readRDS(file.path(processed_dir, inputs[2]))
macro_longterm <- readRDS(file.path(processed_dir, inputs[3]))
taxa_cols <- readRDS(file.path(processed_dir, inputs[4]))

if (!is.character(taxa_cols) || !length(taxa_cols) || anyDuplicated(taxa_cols))
  stop("Taxa column vector is empty, duplicated or not character.")
for (nm in c("macro_metrics", "macro_study_metrics", "macro_longterm")) {
  x <- get(nm)
  missing <- setdiff(c(taxa_cols, "site_id", "working_name", "year",
                       "sample_id", "watercourse", "relative_position",
                       "site_role"), names(x))
  if (length(missing)) stop(nm, " missing fields: ",
                            paste(missing, collapse = ", "))
  if (!nrow(x)) stop(nm, " contains no samples.")
}
if (anyDuplicated(macro_metrics$sample_id) ||
    anyDuplicated(macro_study_metrics$sample_id) ||
    anyDuplicated(macro_longterm$sample_id))
  stop("Duplicate sample_id: inspect before joining EPT metrics.")

# ---- 1. Provisional taxonomic resolution audit --------------------------
# Prefix/underscore classification is a naming heuristic, not verified rank.
taxa_audit <- tibble(taxon = taxa_cols) %>%
  mutate(taxon_level_guess = case_when(
    str_starts(taxon, "o_") ~ "order",
    str_starts(taxon, "f_") ~ "family",
    str_count(taxon, "_") >= 1 &
      !str_starts(taxon, "f_") & !str_starts(taxon, "o_") ~ "genus_or_species",
    TRUE ~ "unclassified"
  ))
taxa_resolution_summary <- taxa_audit %>%
  count(taxon_level_guess, sort = TRUE)
write_csv(taxa_audit, out("step03_taxa_audit.csv"))
write_csv(taxa_resolution_summary, out("step03_taxa_resolution_summary.csv"))

# ---- 2. Provisional EPT candidate list ----------------------------------
# Preserve original v01 keyword list and exclusion for output comparability.
ept_keywords <- c(
  "baetis", "baetidae", "ephemerella", "ephemerellidae",
  "ephemera", "ephemeridae", "ecdyonurus", "heptagenia",
  "heptageniidae", "caenis", "caenidae", "leptophlebiidae",
  "paraleptophlebia", "seratella", "cloeon",
  "plecoptera", "nemouridae", "nemoura", "leuctra", "leuctridae",
  "isoperla", "perlidae", "perlodidae", "chloroperlidae",
  "trichoptera", "hydropsyche", "hydropsychidae",
  "ryacophila", "rhyacophilidae", "limnephilidae",
  "limnephilus", "sericostomatidae", "sericostoma",
  "leptoceridae", "athripsodes", "polycentropodidae",
  "polycentropus", "philopotamidae", "psychomyiidae",
  "glossosomatidae", "goeridae", "lepidosostomatidae",
  "hydroptilidae", "molannidae"
)
candidate_ept <- taxa_cols[str_detect(
  taxa_cols, regex(paste(ept_keywords, collapse = "|"), ignore_case = TRUE)
)]
ept_taxa <- setdiff(candidate_ept, "epoicocladius_ephemerae")
if (!length(ept_taxa)) stop("No provisional EPT taxa identified.")
write_csv(tibble(ept_taxon = ept_taxa),
          out("step03_provisional_ept_taxa.csv"))
write_csv(tibble(
  taxon = candidate_ept,
  included_in_ept = candidate_ept %in% ept_taxa,
  review_status = if_else(candidate_ept %in% ept_taxa,
                          "Provisional: verify taxonomic assignment",
                          "Excluded: known false-positive name match")
), out("step03_ept_candidate_audit.csv"))

# ---- 3. Calculate EPT richness without sample-ID joins ------------------
# Match samples by their complete existing rows; this avoids accidental
# many-to-many sample_id joins and preserves metadata and original row order.
add_ept <- function(x) {
  x %>%
    rowwise() %>%
    mutate(ept_richness = sum(c_across(all_of(ept_taxa)) > 0, na.rm = TRUE)) %>%
    ungroup()
}
macro_metrics_ept <- add_ept(macro_metrics)
macro_study_ept <- add_ept(macro_study_metrics)
macro_longterm_ept <- add_ept(macro_longterm)
stopifnot(nrow(macro_metrics_ept) == nrow(macro_metrics),
          nrow(macro_study_ept) == nrow(macro_study_metrics),
          nrow(macro_longterm_ept) == nrow(macro_longterm),
          all(macro_study_ept$ept_richness <= macro_study_ept$taxon_richness))

# ---- 4. Study-site EPT summaries and plot -------------------------------
study_ept_summary <- macro_study_ept %>%
  group_by(working_name, watercourse, relative_position, site_role) %>%
  summarise(
    first_year = min(year, na.rm = TRUE),
    last_year = max(year, na.rm = TRUE),
    n_samples = n(),
    mean_ept_richness = mean(ept_richness, na.rm = TRUE),
    median_ept_richness = median(ept_richness, na.rm = TRUE),
    min_ept_richness = min(ept_richness, na.rm = TRUE),
    max_ept_richness = max(ept_richness, na.rm = TRUE),
    sd_ept_richness = sd(ept_richness, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(study_ept_summary, out("step03_study_ept_summary.csv"))

p_study_ept_boxplot <- ggplot(
  macro_study_ept, aes(working_name, ept_richness)
) +
  geom_boxplot(outlier.alpha = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  coord_flip() +
  labs(title = "EPT richness across restoration monitoring sites",
       x = "Monitoring site", y = "EPT richness per sample") +
  theme_minimal()
ggsave(out("step03_study_ept_boxplot.png"),
       p_study_ept_boxplot, width = 8, height = 5, dpi = 300)

# ---- 5. Long-term annual EPT summaries and plot -------------------------
longterm_ept_year <- macro_longterm_ept %>%
  group_by(working_name, watercourse, relative_position, site_role, year) %>%
  summarise(
    n_samples = n(),
    mean_ept_richness = mean(ept_richness, na.rm = TRUE),
    median_ept_richness = median(ept_richness, na.rm = TRUE),
    min_ept_richness = min(ept_richness, na.rm = TRUE),
    max_ept_richness = max(ept_richness, na.rm = TRUE),
    .groups = "drop"
  )
p_longterm_ept <- ggplot(longterm_ept_year,
                         aes(year, median_ept_richness)) +
  geom_point(aes(size = n_samples), alpha = 0.8) +
  geom_line(alpha = 0.6) +
  geom_smooth(method = "loess", se = TRUE) +
  geom_vline(xintercept = 2021, colour = "grey70", linetype = "dashed") +
  facet_wrap(~ working_name, scales = "free_y") +
  labs(title = "Long-term EPT richness at core restoration inference sites",
       subtitle = "Dashed line marks the 2021 implementation year",
       x = "Year", y = "Median EPT richness per sample",
       size = "Samples/year") +
  theme_minimal()
ggsave(out("step03_longterm_ept_through_time.png"),
       p_longterm_ept, width = 10, height = 7, dpi = 300)
write_csv(longterm_ept_year, out("step03_longterm_ept_year.csv"))

# ---- 6. Exploratory non-parametric site comparisons ---------------------
ept_kruskal_study <- kruskal.test(ept_richness ~ working_name,
                                 data = macro_study_ept)
ept_pairwise_study <- pairwise.wilcox.test(
  macro_study_ept$ept_richness, macro_study_ept$working_name,
  p.adjust.method = "BH"
)
ept_kruskal_longterm <- kruskal.test(ept_richness ~ working_name,
                                    data = macro_longterm_ept)
ept_pairwise_longterm <- pairwise.wilcox.test(
  macro_longterm_ept$ept_richness, macro_longterm_ept$working_name,
  p.adjust.method = "BH"
)
ept_test_summary <- tibble(
  comparison = c("Six study sites", "Three long-term sites"),
  statistic = c(unname(ept_kruskal_study$statistic),
                unname(ept_kruskal_longterm$statistic)),
  df = c(unname(ept_kruskal_study$parameter),
         unname(ept_kruskal_longterm$parameter)),
  p_value = c(ept_kruskal_study$p.value, ept_kruskal_longterm$p.value)
)
pairwise_table <- function(x, label) {
  as.data.frame(as.table(x$p.value), stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    rename(site_1 = Var1, site_2 = Var2, p_adjusted_BH = Freq) %>%
    filter(!is.na(p_adjusted_BH)) %>%
    mutate(comparison = label, .before = 1)
}
ept_pairwise_summary <- bind_rows(
  pairwise_table(ept_pairwise_study, "Six study sites"),
  pairwise_table(ept_pairwise_longterm, "Three long-term sites")
)
write_csv(ept_test_summary, out("step03_ept_kruskal_tests.csv"))
write_csv(ept_pairwise_summary, out("step03_ept_pairwise_BH.csv"))

# ---- 7. Save processed datasets for downstream analyses -----------------
saveRDS(macro_metrics_ept,
        file.path(processed_dir, "macro_metrics_ept_v01.rds"))
saveRDS(macro_study_ept,
        file.path(processed_dir, "macro_study_ept_v01.rds"))
saveRDS(macro_longterm_ept,
        file.path(processed_dir, "macro_longterm_ept_v01.rds"))

# ---- 8. End-of-script QC ------------------------------------------------
message("Script 04 complete: macro_analysis_step03_v02.R")
message("Provisional EPT columns: ", length(ept_taxa),
        "; study samples: ", nrow(macro_study_ept),
        "; long-term samples: ", nrow(macro_longterm_ept))
message("Six-site EPT Kruskal-Wallis p: ",
        signif(ept_kruskal_study$p.value, 4))
message("Outputs: ", output_dir)
message("Review provisional EPT taxa before ecological interpretation.")
