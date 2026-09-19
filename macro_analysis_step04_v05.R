# ============================================================================
# River Mease | Historical macroinvertebrates
# Script 05: BMWP and ASPT scoring | macro_analysis_step04_v05.R
# ============================================================================
# PURPOSE: Score the six historical study sites using the completed, manually
# reviewed 323-taxon BMWP lookup. Preserve original historical RDS/CSV outputs.
# RUN ORDER: after macro_analysis_step01_v02.R; run from EA_data project root.
# INPUTS: data/processed/macro_study_metrics_v01.rds
#         data/processed/macro_taxa_cols_v01.rds
#         data/processed/bmwp_taxon_lookup_completed_v01.csv
# OUTPUTS: original historical outputs/tables and data/processed filenames.
# SCORING: presence/absence; each distinct scored BMWP group counted once per
# sample, regardless of number of constituent source taxa. Unmapped taxa do not
# score. BMWP=0 and ASPT=NA when no scoring groups occur.
# LIMITATIONS: lookup taxonomy and score assignments are taken as supplied;
# mixed historical taxonomic resolution and sampling effort limit comparison.
# The original script's 2026 GIS appendage depends on external objects and is
# intentionally excluded; it belongs in the separate contemporary GIS workflow.
# STATUS: revised from macro_analysis_step04_v04.R; run and compare against
# submitted results before considering this version validated.
# ============================================================================

library(tidyverse)
processed_dir <- 'data/processed'
output_dir <- 'outputs/tables'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(x) file.path(output_dir, x)
paths <- file.path(processed_dir, c('macro_study_metrics_v01.rds',
                                   'macro_taxa_cols_v01.rds',
                                   'bmwp_taxon_lookup_completed_v01.csv'))
if (any(!file.exists(paths))) stop('Missing input(s): ', paste(paths[!file.exists(paths)], collapse=', '))
macro_study_metrics <- readRDS(paths[1])
taxa_cols <- readRDS(paths[2])
# Ignore the source CSV's unnamed ninth column; preserve all eight scoring fields.
bmwp_taxon_lookup <- read_csv(paths[3], show_col_types=FALSE, name_repair='unique',
                             col_types=cols(.default=col_character())) %>%
  select(dataset_taxon, resolved_taxon, bmwp_group, bmwp_score, mapping_notes) %>%
  mutate(bmwp_score=as.numeric(bmwp_score))

# ---- 1. Dataset and lookup quality checks --------------------------------
stopifnot(nrow(macro_study_metrics)==147,
          n_distinct(macro_study_metrics$site_id)==6,
          length(taxa_cols)>0, !anyDuplicated(taxa_cols),
          all(taxa_cols %in% names(macro_study_metrics)),
          nrow(bmwp_taxon_lookup)==323,
          !anyDuplicated(bmwp_taxon_lookup$dataset_taxon),
          !anyNA(bmwp_taxon_lookup$dataset_taxon))
macro_bmwp <- macro_study_metrics
non_taxa_cols <- c('richness_no_species','total_abundance')
candidate_bmwp_taxa_cols <- setdiff(intersect(taxa_cols, names(macro_bmwp)), non_taxa_cols)
# Taxa are numeric recorded abundances in the processed source; NA = no record.
taxon_occurrence <- macro_bmwp %>%
  summarise(across(all_of(candidate_bmwp_taxa_cols),
                   ~sum(replace_na(as.numeric(.x),0)>0))) %>%
  pivot_longer(everything(), names_to='dataset_taxon', values_to='n_samples')
bmwp_taxa_cols <- taxon_occurrence %>% filter(n_samples>0) %>% pull(dataset_taxon)
missing_lookup <- setdiff(bmwp_taxa_cols, bmwp_taxon_lookup$dataset_taxon)
extra_lookup <- setdiff(bmwp_taxon_lookup$dataset_taxon, bmwp_taxa_cols)
if (length(bmwp_taxa_cols)!=323 || length(missing_lookup) || length(extra_lookup))
  stop('Observed taxa/lookup mismatch. Observed=',length(bmwp_taxa_cols),
       '; missing=',paste(missing_lookup,collapse=', '),
       '; extra=',paste(extra_lookup,collapse=', '))
if (any(!is.na(bmwp_taxon_lookup$bmwp_group) & is.na(bmwp_taxon_lookup$bmwp_score)))
  stop('Mapped BMWP groups missing scores.')
if (any(!is.na(bmwp_taxon_lookup$bmwp_group) &
        (is.na(bmwp_taxon_lookup$bmwp_score) | bmwp_taxon_lookup$bmwp_score<=0)))
  stop('Invalid mapped BMWP score.')
score_conflicts <- bmwp_taxon_lookup %>% filter(!is.na(bmwp_group)) %>%
  distinct(bmwp_group,bmwp_score) %>% count(bmwp_group) %>% filter(n>1)
if (nrow(score_conflicts)) stop('Conflicting scores for BMWP group(s): ',
                                paste(score_conflicts$bmwp_group,collapse=', '))

# ---- 2. Occurrences and operational BMWP groups --------------------------
macro_bmwp_indexed <- macro_bmwp %>% mutate(bmwp_sample_id=row_number())
macro_bmwp_long <- macro_bmwp_indexed %>%
  select(bmwp_sample_id,all_of(bmwp_taxa_cols)) %>%
  pivot_longer(all_of(bmwp_taxa_cols),names_to='dataset_taxon',
               values_to='recorded_abundance') %>%
  mutate(recorded_abundance=replace_na(as.numeric(recorded_abundance),0)) %>%
  filter(recorded_abundance>0) %>%
  left_join(bmwp_taxon_lookup,by='dataset_taxon',relationship='many-to-one')
bmwp_unscored_occurrences <- macro_bmwp_long %>% filter(is.na(bmwp_group)) %>%
  count(dataset_taxon,resolved_taxon,mapping_notes,sort=TRUE,
        name='n_sample_occurrences')
unexplained <- macro_bmwp_long %>%
  filter(is.na(bmwp_group),is.na(mapping_notes) | trimws(mapping_notes)=='') %>%
  distinct(dataset_taxon)
if (nrow(unexplained)) stop('Unmapped taxa without explanatory notes: ',
                           paste(unexplained$dataset_taxon,collapse=', '))
bmwp_group_presence <- macro_bmwp_long %>%
  filter(!is.na(bmwp_group),!is.na(bmwp_score)) %>%
  distinct(bmwp_sample_id,bmwp_group,bmwp_score)
bmwp_group_presence_metadata <- bmwp_group_presence %>%
  left_join(macro_bmwp_indexed %>%
              select(bmwp_sample_id,site_id,working_name,watercourse,
                     relative_position,site_role,sample_date,year),
            by='bmwp_sample_id',relationship='many-to-one') %>%
  select(bmwp_sample_id,site_id,working_name,watercourse,relative_position,
         site_role,sample_date,year,bmwp_group,bmwp_score) %>%
  arrange(working_name,sample_date,bmwp_group)

# ---- 3. Score once per group, per sample ---------------------------------
bmwp_sample_scores <- bmwp_group_presence %>% group_by(bmwp_sample_id) %>%
  summarise(bmwp=sum(bmwp_score),n_bmwp_groups=n_distinct(bmwp_group),
            aspt=bmwp/n_bmwp_groups,.groups='drop')
macro_bmwp_scored <- macro_bmwp_indexed %>%
  left_join(bmwp_sample_scores,by='bmwp_sample_id',relationship='one-to-one') %>%
  mutate(bmwp=replace_na(bmwp,0),
         n_bmwp_groups=replace_na(n_bmwp_groups,0L)) %>%
  select(-bmwp_sample_id)
stopifnot(nrow(macro_bmwp_scored)==147,
          n_distinct(macro_bmwp_scored$site_id)==6,
          all(macro_bmwp_scored$bmwp>=0),
          all(is.na(macro_bmwp_scored$aspt)==(macro_bmwp_scored$n_bmwp_groups==0)))

# ---- 4. Historical exports (original names) ------------------------------
saveRDS(macro_bmwp_scored,file.path(processed_dir,'macro_bmwp_scored_v01.rds'))
saveRDS(macro_bmwp_long,file.path(processed_dir,'macro_bmwp_long_v01.rds'))
saveRDS(bmwp_group_presence,file.path(processed_dir,'bmwp_group_presence_v01.rds'))
saveRDS(bmwp_group_presence_metadata,
        file.path(processed_dir,'bmwp_group_presence_metadata_v01.rds'))
write_csv(macro_bmwp_scored,out('macro_bmwp_scored_v01.csv'),na='')
write_csv(bmwp_unscored_occurrences,out('bmwp_unscored_occurrences_v01.csv'),na='')
write_csv(taxon_occurrence,out('bmwp_taxon_occurrence_v01.csv'),na='')
write_csv(bmwp_group_presence_metadata,out('bmwp_group_presence_metadata_v01.csv'),na='')
historic_bmwp_group_gis <- macro_bmwp_scored %>%
  group_by(working_name,watercourse,relative_position,site_role) %>%
  summarise(n_samples=n(),median_bmwp_groups=median(n_bmwp_groups),
            mean_bmwp_groups=mean(n_bmwp_groups),
            min_bmwp_groups=min(n_bmwp_groups),
            max_bmwp_groups=max(n_bmwp_groups),.groups='drop')
write_csv(historic_bmwp_group_gis,out('historic_bmwp_group_gis_v01.csv'))

# Additional compact audit outputs for comparison with thesis.
site_bmwp_summary <- macro_bmwp_scored %>%
  group_by(working_name,site_id,site_role) %>%
  summarise(n_samples=n(),mean_bmwp=mean(bmwp),median_bmwp=median(bmwp),
            sd_bmwp=sd(bmwp),mean_aspt=mean(aspt,na.rm=TRUE),
            median_aspt=median(aspt,na.rm=TRUE),sd_aspt=sd(aspt,na.rm=TRUE),
            mean_n_bmwp_groups=mean(n_bmwp_groups),.groups='drop')
write_csv(site_bmwp_summary,out('historic_bmwp_site_summary_v01.csv'))
cat('\n=== Historical BMWP/ASPT scoring complete ===\n')
cat('Samples:',nrow(macro_bmwp_scored),'| Observed source taxa:',
    length(bmwp_taxa_cols),'| Unmapped taxa:',
    sum(is.na(bmwp_taxon_lookup$bmwp_group)),'\n')
print(site_bmwp_summary,n=Inf,width=Inf)
cat('Historical outputs written to',output_dir,'and',processed_dir,'\n')
cat('2026 GIS exports are handled separately; no external 2026 objects required.\n')
