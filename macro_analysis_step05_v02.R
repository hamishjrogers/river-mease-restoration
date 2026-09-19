# ============================================================================
# River Mease | Historical macroinvertebrates
# Script 06: Historical BMWP and ASPT summaries, site tests and figures
# File: macro_analysis_step05_v02.R
# ============================================================================
# RUN AFTER: macro_analysis_step04_v05.R; working directory = EA_data root.
# INPUT: data/processed/macro_bmwp_scored_v01.rds (historical samples only).
# OUTPUT: original step05 CSV and PNG filenames, plus QC summary.
# METHOD: annual and site-level descriptive statistics; exploratory Kruskal-
# Wallis and pairwise Wilcoxon tests (BH adjustment; exact = FALSE). Unequal
# monitoring effort and repeated samples limit independent-sample inference.
# The 2021 line is a chronological marker, NOT a causal restoration estimate.
# Original step05 script ends partway through its publication-figure section;
# this version supplies both long-term figures with gaps between sampled years.
# STATUS: revised code; run locally and compare with submitted thesis figures.
# ============================================================================

library(tidyverse)
input_file <- 'data/processed/macro_bmwp_scored_v01.rds'
if (!file.exists(input_file)) stop('Missing ', input_file, '; run Script 05 first.')
dir.create('outputs/tables', recursive = TRUE, showWarnings = FALSE)
dir.create('outputs/figures', recursive = TRUE, showWarnings = FALSE)
table_out <- function(x) file.path('outputs/tables', x)
figure_out <- function(x) file.path('outputs/figures', x)

site_order <- c('Upper Gilwiskaw', 'B-road Gilwiskaw', 'Upper Mease',
                'Restoration Reach Mease', 'Birds Hill', 'Lower Mease')
longterm_site_order <- c('Upper Gilwiskaw', 'Restoration Reach Mease',
                         'Lower Mease')
site_letters <- c('Upper Gilwiskaw'='A', 'B-road Gilwiskaw'='B',
                  'Upper Mease'='C', 'Restoration Reach Mease'='D',
                  'Birds Hill'='E', 'Lower Mease'='F')
macro_bmwp_scored <- readRDS(input_file)
required <- c('working_name', 'site_id', 'site_role', 'watercourse',
              'relative_position', 'year', 'bmwp', 'aspt', 'n_bmwp_groups')
missing <- setdiff(required, names(macro_bmwp_scored))
if (length(missing)) stop('Missing fields: ', paste(missing, collapse=', '))
if (nrow(macro_bmwp_scored) != 147L ||
    n_distinct(macro_bmwp_scored$site_id) != 6L ||
    anyNA(macro_bmwp_scored$year) ||
    anyNA(macro_bmwp_scored$working_name) ||
    !setequal(unique(as.character(macro_bmwp_scored$working_name)), site_order))
  stop('Historical sample count, site identities or year/label QC failed.')
if (anyNA(macro_bmwp_scored$bmwp) || anyNA(macro_bmwp_scored$aspt))
  stop('Missing BMWP or ASPT: inspect scoring output before testing.')
macro_bmwp_scored <- macro_bmwp_scored %>%
  mutate(working_name = factor(as.character(working_name), levels=site_order))

# ---- 1. Annual descriptive summaries -----------------------------------
annual_bmwp_summary <- macro_bmwp_scored %>%
  group_by(working_name, watercourse, relative_position, site_role, year) %>%
  summarise(n_samples=sum(!is.na(bmwp)), mean_bmwp=mean(bmwp,na.rm=TRUE),
            median_bmwp=median(bmwp,na.rm=TRUE), sd_bmwp=sd(bmwp,na.rm=TRUE),
            min_bmwp=min(bmwp,na.rm=TRUE), max_bmwp=max(bmwp,na.rm=TRUE),
            .groups='drop') %>% arrange(working_name,year)
annual_aspt_summary <- macro_bmwp_scored %>%
  group_by(working_name, watercourse, relative_position, site_role, year) %>%
  summarise(n_samples=sum(!is.na(aspt)), mean_aspt=mean(aspt,na.rm=TRUE),
            median_aspt=median(aspt,na.rm=TRUE), sd_aspt=sd(aspt,na.rm=TRUE),
            min_aspt=min(aspt,na.rm=TRUE), max_aspt=max(aspt,na.rm=TRUE),
            .groups='drop') %>% arrange(working_name,year)
write_csv(annual_bmwp_summary,table_out('annual_bmwp_summary_v01.csv'))
write_csv(annual_aspt_summary,table_out('annual_aspt_summary_v01.csv'))

# ---- 2. Overall site summaries (preserve original column names) --------
bmwp_site_summary <- macro_bmwp_scored %>%
  group_by(working_name, watercourse, relative_position, site_role) %>%
  summarise(n=n(), sample_period=paste0(min(year),'–',max(year)),
            median_bmwp=median(bmwp,na.rm=TRUE),
            mean_bmwp=mean(bmwp,na.rm=TRUE), sd_bmwp=sd(bmwp,na.rm=TRUE),
            range_bmwp=paste0(min(bmwp,na.rm=TRUE),'–',max(bmwp,na.rm=TRUE)),
            .groups='drop') %>% arrange(working_name)
aspt_site_summary <- macro_bmwp_scored %>%
  group_by(working_name, watercourse, relative_position, site_role) %>%
  summarise(n=n(), sample_period=paste0(min(year,na.rm=TRUE),'–',max(year,na.rm=TRUE)),
            median_aspt=median(aspt,na.rm=TRUE),
            mean_aspt=mean(aspt,na.rm=TRUE), sd_aspt=sd(aspt,na.rm=TRUE),
            range_aspt=paste0(round(min(aspt,na.rm=TRUE),2),'–',
                              round(max(aspt,na.rm=TRUE),2)),
            .groups='drop') %>% arrange(working_name)
write_csv(bmwp_site_summary,table_out('bmwp_site_summary_v01.csv'))
write_csv(aspt_site_summary,table_out('aspt_site_summary_v01.csv'))

# ---- 3. Exploratory site comparisons ------------------------------------
bmwp_kw <- kruskal.test(bmwp ~ working_name, data=macro_bmwp_scored)
aspt_kw <- kruskal.test(aspt ~ working_name, data=macro_bmwp_scored)
bmwp_aspt_kw_results <- tibble(
  response=c('BMWP','ASPT'),
  statistic=c(unname(bmwp_kw$statistic),unname(aspt_kw$statistic)),
  degrees_freedom=c(unname(bmwp_kw$parameter),unname(aspt_kw$parameter)),
  p_value=c(bmwp_kw$p.value,aspt_kw$p.value))
make_pairwise <- function(metric) {
  test <- pairwise.wilcox.test(macro_bmwp_scored[[metric]],
                               macro_bmwp_scored$working_name,
                               p.adjust.method='BH',exact=FALSE)
  as.data.frame(as.table(test$p.value),stringsAsFactors=FALSE) %>%
    as_tibble() %>% rename(site_1=Var1,site_2=Var2,adjusted_p_value=Freq) %>%
    filter(!is.na(adjusted_p_value)) %>% arrange(adjusted_p_value)
}
bmwp_pairwise_table <- make_pairwise('bmwp')
aspt_pairwise_table <- make_pairwise('aspt')
bmwp_pairwise_significant <- bmwp_pairwise_table %>%
  filter(adjusted_p_value < 0.05)
aspt_pairwise_significant <- aspt_pairwise_table %>%
  filter(adjusted_p_value < 0.05)
write_csv(bmwp_aspt_kw_results,table_out('bmwp_aspt_kruskal_results_v01.csv'))
write_csv(bmwp_pairwise_table,table_out('bmwp_pairwise_wilcoxon_v01.csv'))
write_csv(aspt_pairwise_table,table_out('aspt_pairwise_wilcoxon_v01.csv'))

# ---- 4. Six-site boxplots and annual time series ------------------------
plot_order <- rev(site_order)
boxplot_metric <- function(metric,label) {
  macro_bmwp_scored %>%
    mutate(working_name=factor(working_name,levels=plot_order)) %>%
    ggplot(aes(y=working_name,x=.data[[metric]])) +
    geom_boxplot(linewidth=0.8,outlier.shape=NA) +
    geom_jitter(height=0.15,alpha=0.35,size=2) +
    scale_y_discrete(labels=site_letters) +
    scale_x_continuous(expand=expansion(mult=c(0.02,0.04))) +
    labs(x=label,y='Site') + theme_bw(base_size=12) +
    theme(panel.grid.minor=element_blank(),
          axis.title=element_text(face='bold'),
          axis.text=element_text(colour='black'))
}
bmwp_boxplot <- boxplot_metric('bmwp','BMWP score')
aspt_boxplot <- boxplot_metric('aspt','ASPT')
plot_annual <- function(x, metric, label, title) {
  ggplot(x,aes(year,.data[[metric]],group=working_name)) +
    geom_line(linewidth=0.6,na.rm=TRUE) + geom_point(size=2,na.rm=TRUE) +
    facet_wrap(~working_name,scales='free_x',ncol=2) +
    geom_vline(xintercept=2021,linetype='dashed',linewidth=0.5) +
    labs(x='Year',y=label,title=title,
         subtitle='Dashed line indicates restoration implementation in 2021') +
    theme_bw() + theme(plot.title=element_text(face='bold'),
                       panel.grid.minor=element_blank())
}
bmwp_time_series <- plot_annual(annual_bmwp_summary,'median_bmwp',
                                 'Annual median BMWP','BMWP through time at the six study sites')
aspt_time_series <- plot_annual(annual_aspt_summary,'median_aspt',
                                 'Annual median ASPT','ASPT through time at the six study sites')
ggsave(figure_out('bmwp_boxplot_by_site_v01.png'),bmwp_boxplot,width=9,height=6,dpi=300)
ggsave(figure_out('aspt_boxplot_by_site_v01.png'),aspt_boxplot,width=9,height=6,dpi=300)
ggsave(figure_out('bmwp_time_series_v01.png'),bmwp_time_series,width=10,height=10,dpi=300)
ggsave(figure_out('aspt_time_series_v01.png'),aspt_time_series,width=10,height=10,dpi=300)

# ---- 5. Long-term publication figures ----------------------------------
# Unlike the simple six-site plots, connect only consecutive sampled years.
prepare_longterm <- function(x) {
  x %>% filter(as.character(working_name) %in% longterm_site_order) %>%
    mutate(working_name=factor(as.character(working_name),
                               levels=longterm_site_order)) %>%
    arrange(working_name,year) %>% group_by(working_name) %>%
    mutate(year_gap=year-lag(year),
           segment_start=if_else(is.na(year_gap)|year_gap>1,1L,0L),
           line_segment=cumsum(segment_start)) %>% ungroup()
}
bmwp_longterm_plot_data <- prepare_longterm(annual_bmwp_summary)
aspt_longterm_plot_data <- prepare_longterm(annual_aspt_summary)
plot_publication <- function(x,metric,label) {
  ggplot(x,aes(year,.data[[metric]])) +
    geom_line(aes(group=interaction(working_name,line_segment)),
              linewidth=0.65,na.rm=TRUE) +
    geom_point(aes(size=n_samples),shape=21,fill='white',stroke=0.7,na.rm=TRUE) +
    geom_vline(xintercept=2021,linetype='dashed',linewidth=0.6) +
    facet_wrap(~working_name,ncol=1,
      labeller=as_labeller(c('Upper Gilwiskaw'='Site A',
                            'Restoration Reach Mease'='Site D',
                            'Lower Mease'='Site F'))) +
    scale_x_continuous(breaks=seq(1985,2025,5),limits=c(1985,2025),
                       expand=expansion(mult=c(0.01,0.02))) +
    scale_size_continuous(name='Samples per year',range=c(2.2,4.2),
                          breaks=c(1,2,3,4)) +
    labs(x='Year',y=label) + theme_bw(base_size=11) +
    theme(strip.background=element_rect(fill='grey95',colour='grey40',linewidth=0.4),
          strip.text=element_text(face='bold',size=10.5),
          panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),
          panel.spacing=grid::unit(0.8,'lines'),
          axis.title=element_text(face='bold'),legend.position='right')
}
bmwp_time_series_publication <- plot_publication(
  bmwp_longterm_plot_data,'median_bmwp','Annual median BMWP score')
aspt_time_series_publication <- plot_publication(
  aspt_longterm_plot_data,'median_aspt','Annual median ASPT')
ggsave(figure_out('bmwp_time_series_publication_v01.png'),
       bmwp_time_series_publication,width=8,height=9,dpi=300)
ggsave(figure_out('aspt_time_series_publication_v01.png'),
       aspt_time_series_publication,width=8,height=9,dpi=300)

# ---- 6. End-of-script checks --------------------------------------------
cat('\n=== Historical BMWP/ASPT site comparisons complete ===\n')
cat('Historical samples:',nrow(macro_bmwp_scored),'| Sites:',
    n_distinct(macro_bmwp_scored$site_id),'\n')
print(bmwp_aspt_kw_results,n=Inf,width=Inf)
cat('\nSignificant BH-adjusted BMWP comparisons:\n')
print(bmwp_pairwise_significant,n=Inf,width=Inf)
cat('\nSignificant BH-adjusted ASPT comparisons:\n')
print(aspt_pairwise_significant,n=Inf,width=Inf)
cat('\nOutputs: outputs/tables and outputs/figures\n')
