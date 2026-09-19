# ============================================================================
# River Mease | Contemporary macroinvertebrates
# Script 08: Descriptive BMWP and ASPT at six 2026 field sites
# File: macro_2026BMWPanalysis_v02.R
# ============================================================================
# INPUT: data/processed/macro_bmwp_scored_combined_2026_v01.rds (Script 07)
# RUN FROM: /home/hamish/Documents/NTU/research_project/data_analysis/data/EA_data
# PURPOSE: Describe the single 19 July 2026 sample at each of six field sites.
# No inferential comparisons: there is no within-site replication.
# Preserve existing descriptive output filenames; community-dissimilarity
# appendix tables B3.8 and B3.9 belong to the later community analysis.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

input_path <- 'data/processed/macro_bmwp_scored_combined_2026_v01.rds'
if (!file.exists(input_path)) stop('Missing Script 07 output: ', input_path)
dir.create('outputs/tables', recursive = TRUE, showWarnings = FALSE)
dir.create('outputs/figures', recursive = TRUE, showWarnings = FALSE)

macro_bmwp_scored_combined <- readRDS(input_path)
required <- c('field_site','working_name','year','data_source','sample_date',
              'bmwp','aspt','n_bmwp_groups')
missing <- setdiff(required, names(macro_bmwp_scored_combined))
if (length(missing)) stop('Missing required columns: ', paste(missing, collapse = ', '))
if (nrow(macro_bmwp_scored_combined) != 153L ||
    sum(macro_bmwp_scored_combined$data_source == 'Historic secondary data', na.rm = TRUE) != 147L)
  stop('Expected 153 combined rows including 147 historical samples.')

field_site_order <- c('Downstream Mease','Confluence','Wetland Mease',
                      'Upstream Mease','Wetland Gilwiskaw','Upstream Gilwiskaw')
site_number_lookup <- c('Upstream Gilwiskaw'=1L,'Wetland Gilwiskaw'=2L,
                        'Upstream Mease'=3L,'Wetland Mease'=4L,
                        'Confluence'=5L,'Downstream Mease'=6L)

# Use the original field_site field, not working_name: the latter intentionally
# maps two contemporary locations to nearby historical monitoring site labels.
bmwp_aspt_2026 <- macro_bmwp_scored_combined %>%
  filter(data_source == '2026 field survey', year == 2026) %>%
  mutate(field_site = as.character(field_site),
         spatial_group = case_when(
           field_site %in% c('Upstream Gilwiskaw','Upstream Mease') ~ 'Upstream',
           field_site %in% c('Wetland Gilwiskaw','Wetland Mease') ~ 'Wetland',
           field_site %in% c('Confluence','Downstream Mease') ~ 'Downstream',
           TRUE ~ NA_character_),
         reach_type = case_when(
           field_site == 'Upstream Gilwiskaw' ~ 'Upstream tributary',
           field_site == 'Wetland Gilwiskaw' ~ 'Restored wetland tributary',
           field_site == 'Upstream Mease' ~ 'Upstream main channel',
           field_site == 'Wetland Mease' ~ 'Restored wetland main channel',
           field_site == 'Confluence' ~ 'River confluence',
           field_site == 'Downstream Mease' ~ 'Downstream main channel',
           TRUE ~ NA_character_),
         watercourse_2026 = case_when(
           field_site %in% c('Upstream Gilwiskaw','Wetland Gilwiskaw') ~ 'Gilwiskaw',
           field_site %in% c('Upstream Mease','Wetland Mease','Downstream Mease') ~ 'Mease',
           field_site == 'Confluence' ~ 'Confluence',
           TRUE ~ NA_character_))

if (nrow(bmwp_aspt_2026) != 6L ||
    !setequal(bmwp_aspt_2026$field_site, field_site_order) ||
    anyDuplicated(bmwp_aspt_2026$field_site) ||
    anyNA(bmwp_aspt_2026[,c('bmwp','aspt','n_bmwp_groups',
                            'spatial_group','reach_type','watercourse_2026')]) ||
    any(bmwp_aspt_2026$n_bmwp_groups <= 0) ||
    any(abs(bmwp_aspt_2026$bmwp / bmwp_aspt_2026$n_bmwp_groups -
            bmwp_aspt_2026$aspt) > 1e-8))
  stop('2026 site identities, scores, or BMWP/ASPT arithmetic failed QC.')
if (any(as.Date(bmwp_aspt_2026$sample_date) != as.Date('2026-07-19')))
  stop('2026 macroinvertebrate survey date differs from confirmed field log: 2026-07-19.')

bmwp_aspt_2026_summary <- bmwp_aspt_2026 %>%
  transmute(field_site = factor(field_site, levels = field_site_order),
            spatial_group = factor(spatial_group, levels = c('Upstream','Wetland','Downstream')),
            reach_type, watercourse = watercourse_2026,
            bmwp = as.numeric(bmwp), aspt = as.numeric(aspt),
            n_bmwp_groups = as.integer(n_bmwp_groups)) %>%
  arrange(field_site)

# Rankings are descriptive ordering of single samples, not ecological ratings.
bmwp_ranking_2026 <- bmwp_aspt_2026_summary %>%
  arrange(desc(bmwp)) %>% mutate(bmwp_rank = row_number()) %>%
  select(bmwp_rank, field_site, bmwp)
aspt_ranking_2026 <- bmwp_aspt_2026_summary %>%
  arrange(desc(aspt)) %>% mutate(aspt_rank = row_number()) %>%
  select(aspt_rank, field_site, aspt)

bmwp_aspt_2026_long <- bmwp_aspt_2026_summary %>%
  select(field_site, spatial_group, bmwp, aspt) %>%
  pivot_longer(c(bmwp, aspt), names_to = 'metric', values_to = 'value') %>%
  mutate(metric = factor(recode(metric, bmwp='BMWP score', aspt='ASPT'),
                         levels = c('BMWP score','ASPT')))

shape_values <- c('Upstream'=16, 'Wetland'=17, 'Downstream'=15)
plot_theme <- theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(),
        axis.text=element_text(colour='black'), legend.position='bottom')

bmwp_aspt_2026_plot <- ggplot(bmwp_aspt_2026_long,
                              aes(x=value,y=field_site,shape=spatial_group)) +
  geom_segment(aes(x=0,xend=value,yend=field_site),
               linewidth=0.6,alpha=0.55,show.legend=FALSE) +
  geom_point(size=3.5,stroke=0.8) +
  facet_wrap(~metric,scales='free_x',nrow=1) +
  scale_shape_manual(values=shape_values,name='Spatial setting') +
  scale_x_continuous(expand=expansion(mult=c(0,0.08))) +
  labs(x=NULL,y='2026 field-survey site',
       title='BMWP and ASPT across the 2026 survey sites',
       subtitle='One macroinvertebrate sample per site, 19 July 2026') + plot_theme

bmwp_2026_plot <- ggplot(bmwp_aspt_2026_summary,
                         aes(x=bmwp,y=field_site,shape=spatial_group)) +
  geom_segment(aes(x=0,xend=bmwp,yend=field_site),
               linewidth=0.7,alpha=0.55,show.legend=FALSE) +
  geom_point(size=3.8) + scale_shape_manual(values=shape_values,name='Spatial setting') +
  labs(x='BMWP score',y='2026 field-survey site') + plot_theme
aspt_2026_plot <- ggplot(bmwp_aspt_2026_summary,
                         aes(x=aspt,y=field_site,shape=spatial_group)) +
  geom_segment(aes(x=0,xend=aspt,yend=field_site),
               linewidth=0.7,alpha=0.55,show.legend=FALSE) +
  geom_point(size=3.8) + scale_shape_manual(values=shape_values,name='Spatial setting') +
  labs(x='ASPT',y='2026 field-survey site') + plot_theme

write_csv(bmwp_aspt_2026_summary,'outputs/tables/bmwp_aspt_2026_site_summary_v01.csv')
write_csv(bmwp_ranking_2026,'outputs/tables/bmwp_2026_site_ranking_v01.csv')
write_csv(aspt_ranking_2026,'outputs/tables/aspt_2026_site_ranking_v01.csv')
ggsave('outputs/figures/bmwp_aspt_2026_by_site_v01.png',
       bmwp_aspt_2026_plot,width=10,height=6,dpi=300,bg='white')
ggsave('outputs/figures/bmwp_2026_by_site_v01.png',
       bmwp_2026_plot,width=8,height=5.5,dpi=300,bg='white')
ggsave('outputs/figures/aspt_2026_by_site_v01.png',
       aspt_2026_plot,width=8,height=5.5,dpi=300,bg='white')

appendix_B3_7 <- bmwp_aspt_2026_summary %>%
  mutate(Site = unname(site_number_lookup[as.character(field_site)])) %>%
  transmute(Site, `Site name`=as.character(field_site),
            `Spatial setting`=as.character(spatial_group),
            Watercourse=watercourse,BMWP=bmwp,ASPT=round(aspt,2),
            `BMWP group richness`=n_bmwp_groups) %>% arrange(Site)
write_csv(appendix_B3_7,'outputs/appendix_B3_7_contemporary_metrics.csv')

cat('\n=== 2026 descriptive BMWP/ASPT analysis complete ===\n')
cat('Survey date: 2026-07-19 | Samples:',nrow(bmwp_aspt_2026_summary),'\n')
print(appendix_B3_7,n=Inf,width=Inf)
cat('Saved three summary/ranking tables, three figures and appendix B3.7.\n')
cat('Appendices B3.8 and B3.9 require community-dissimilarity objects and are deferred.\n')
