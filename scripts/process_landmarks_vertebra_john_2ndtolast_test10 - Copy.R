#### Initialize ####
if (
  interactive() &&
    requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()
) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# Install once, if needed:
# install.packages("geomorph")

library(geomorph)
library(tidyverse)
library(janitor)
library(patchwork)

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)




source("gm_functions.R")
getwd()
#### Get Data ####
# Path to your file (adjust if needed)
# tps_file <- "vertebra_2nd-to-last.TPS"


tps_files <-
  list.files(
    path = "../xray_posterior_john_mcphail_non_dvc",
    pattern = "(?i)posterior.*\\.TPS$",
    full.names = TRUE
  ) %>%
  normalizePath()



#### 1. Set up and read the TPS file ####
# Your tpsDig file uses IMAGE= lines (no ID=), so:
# Y <- readland.tps(tps_file,
#                   specID    = "imageID",  # use the IMAGE= names as specimen names
#                   negNA     = TRUE,       # treat negative coords as missing if present
#                   readcurves = FALSE)     # no CURVES= blocks in your file
# dim(Y)

Y_list <- 
  tps_files %>%
  purrr::map(
    function(file) {
      Yi <- readland.tps(
        file,
        specID    = "imageID",
        negNA     = TRUE,
        readcurves = FALSE
      )
      
      # prepend TPS filename to specimen names so you keep file info
      dimnames(Yi)[[3]] <- paste0(basename(file), ":", dimnames(Yi)[[3]])
      
      Yi
    }
  )

Y <- abind::abind(Y_list, along = 3)

dim(Y)

#### 2. Generalized Procrustes Analysis (GPA) ####
Y.gpa <- gpagen(Y,
                print.progress = TRUE,  # echo iteration progress
                ProcD = FALSE)          # just do GPA, no ANOVA here

str(Y.gpa)

#### 3. Outline: Landmark Connections ####

# landmark connections as edges (updated for this image)
outline_edges <-
  tribble(
    ~from, ~to,
    # top loop: 1 -> 5 -> 6 -> 4 -> 3 -> 2 -> 1
    1,    5,
    5,    6,
    6,    4,
    4,    3,
    3,    2,
    2,    1,
    
    # bottom loop: 9 -> 12 -> 7 -> 8 -> 13 -> 11 -> 10 -> 9
    9,    12,
    12,   7,
    7,    8,
    8,    13,
    13,   11,
    11,   10,
    10,   9
  )


#### 4. Turn GPA coordinates into a tidy tibble ####

coords_tbl <- gm_step4_tidy_coords(Y.gpa$coords)

coords_tbl

#### 5. Mean shape + all specimens with ggplot ####

mean_shape_results <- gm_step5_mean_shape_plot(coords_tbl, outline_edges)
mean_shape_tbl <- mean_shape_results$mean_shape_tbl
outline_mean_tbl <- mean_shape_results$outline_mean_tbl
mean_shape_plot <- mean_shape_results$plot +
  scale_x_reverse() +
  scale_y_reverse()

mean_shape_plot



#### 6. Mean shape by era x site with ggplot ####

mean_shape_era_site_results <- gm_step6_mean_shape_by_group(coords_tbl, outline_edges)
mean_shape_tbl_era_site <- mean_shape_era_site_results$mean_shape_tbl
outline_mean_tbl_era_site <- mean_shape_era_site_results$outline_mean_tbl
mean_shape_plot_era_site <- mean_shape_era_site_results$plot

#### 7. PCA with gm.prcomp + tidyverse ####

pca_results <- gm_step7_pca(Y.gpa$coords, gm_parse_specimen)
pca <- pca_results$pca
scores_tbl <- pca_results$scores_tbl
plot_pca <- pca_results$plot_pca
plot_pca_34 <- pca_results$plot_pca_34
plot_pca
plot_pca_34
#### 8. Make Shapes at Ends of PC Axes####

shape_pc_plots <- gm_step8_pc_end_shapes(
  pca = pca,
  outline_edges = outline_edges,
  mean_shape_tbl = mean_shape_tbl,
  outline_mean_tbl = outline_mean_tbl
)

#### 9. Corner shapes in PC1–PC2 space ####

shape_corner_plots <- gm_step9_corner_shapes(
  pca = pca,
  consensus_mat = Y.gpa$consensus,
  outline_edges = outline_edges,
  mean_shape_tbl = mean_shape_tbl,
  outline_mean_tbl = outline_mean_tbl
)

#### 10. Arrange shapes around the PCA plot ####

pca_layouts <- gm_step10_arrange_pca_shapes(
  plot_pca = plot_pca,
  shapes_pc = shape_pc_plots,
  corner_shapes = shape_corner_plots
)

#### 11.1 Allometric effect on shape: Build a geomorph data frame  ####

geomorph_inputs <- gm_step11_1_geomorph_df(
  coords_array = Y.gpa$coords,
  Csize = Y.gpa$Csize,
  coords_tbl = coords_tbl
)

gdf <- geomorph_inputs$gdf
specimen_order <- geomorph_inputs$specimen_order
spec_meta <- geomorph_inputs$spec_meta

#### 11.2. Test for overall allometry: shape ~ log(size) ####

fit_allom <- gm_step11_2_test_allometry(gdf)

#### 11.3. Does allometry differ by era or site? ####

fit_allom_group <- gm_step11_3_test_allometry_by_group(gdf)

#### 11.4. Visualizing the allometric pattern ####

allometry_outputs <- gm_step11_4_allometry_plot(fit_allom, gdf)
reg_scores <- allometry_outputs$reg_scores
allom_df <- allometry_outputs$allom_df
allometry_plot <- allometry_outputs$plot
allometry_plot

#### 11.5 size-corrected shapes (at mean size) ####

shape_size_adjusted <- gm_step11_5_size_corrected_shapes(fit_allom, Y.gpa$coords)

#### 11.6 Test era and site on size-corrected shapes ####

residual_tests <- gm_step11_6_test_residual_shapes(shape_size_adjusted, spec_meta)
gdf_resid <- residual_tests$gdf_resid
fit_resid_era_site <- residual_tests$fit_resid_era_site

#### Repeat steps 4–10 using size-corrected shapes ####

coords_tbl_resid <- gm_step4_tidy_coords(shape_size_adjusted)

mean_shape_results_resid <- gm_step5_mean_shape_plot(
  coords_tbl_resid,
  outline_edges,
  title = "Mean shape (size-corrected, red) with all aligned specimens"
)
mean_shape_tbl_resid <- mean_shape_results_resid$mean_shape_tbl
outline_mean_tbl_resid <- mean_shape_results_resid$outline_mean_tbl
mean_shape_plot_resid <- mean_shape_results_resid$plot
mean_shape_plot_resid

mean_shape_era_site_results_resid <- gm_step6_mean_shape_by_group(coords_tbl_resid, outline_edges)
mean_shape_tbl_era_site_resid <- mean_shape_era_site_results_resid$mean_shape_tbl
outline_mean_tbl_era_site_resid <- mean_shape_era_site_results_resid$outline_mean_tbl
mean_shape_plot_era_site_resid <- mean_shape_era_site_results_resid$plot
mean_shape_plot_era_site_resid

pca_results_resid <- gm_step7_pca(
  shape_size_adjusted,
  gm_parse_specimen,
  title = "PCA of size corrected vertebra shape (PC1 vs PC2)"
)
pca_resid <- pca_results_resid$pca
scores_tbl_resid <- pca_results_resid$scores_tbl
plot_pca_resid <- pca_results_resid$plot_pca
plot_pca_resid_34 <- pca_results_resid$plot_pca_34

shape_pc_plots_resid <- gm_step8_pc_end_shapes(
  pca = pca_resid,
  outline_edges = outline_edges,
  mean_shape_tbl = mean_shape_tbl_resid,
  outline_mean_tbl = outline_mean_tbl_resid
)

consensus_resid <- geomorph::mshape(shape_size_adjusted)

shape_corner_plots_resid <- gm_step9_corner_shapes(
  pca = pca_resid,
  consensus_mat = consensus_resid,
  outline_edges = outline_edges,
  mean_shape_tbl = mean_shape_tbl_resid,
  outline_mean_tbl = outline_mean_tbl_resid
)

pca_layouts_resid <- gm_step10_arrange_pca_shapes(
  plot_pca = plot_pca_resid,
  shapes_pc = shape_pc_plots_resid,
  corner_shapes = shape_corner_plots_resid
)
pca_layouts_resid
