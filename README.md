R scripts supporting the water-quality component of an MSc Biodiversity Conservation dissertation investigating ecological responses to floodplain restoration in the River Mease catchment.

Script execution order

Run the scripts in the following order, from the project root directory:

wetland_restoration_cleaning_v03.R – Cleans the 2026 field monitoring dataset.
trt_2021_diagnostics_v04.R – Examines uncertainty in the reporting units of the 2021 Environment Agency phosphorus data and records the conversion assumption.
trt_cleaning_v06.R – Harmonises historical Environment Agency data with the cleaned 2026 monitoring dataset.
analysis_catchment_context_v03.R – Produces historical catchment-scale phosphorus summaries, statistical analyses and figures.
restoration_analysis_v05.R – Analyses contemporary water-quality conditions, including spatial and temporal variation and SAC phosphorus target exceedances.
trt_sensitivity_analysis_v04.R – Evaluates the sensitivity of historical phosphorus results to alternative reporting-unit assumptions.
Data and reproducibility

The scripts use relative file paths and expect the required input datasets in a data/ directory. Generated analytical outputs and figures are saved to the locations specified within each script.

The underlying datasets are not included in this repository. Reproducing the analyses requires access to the original input data and the R packages loaded by the scripts.
