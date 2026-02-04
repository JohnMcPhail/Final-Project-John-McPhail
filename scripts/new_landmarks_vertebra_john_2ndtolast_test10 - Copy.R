#### Initialize ####
if (
  interactive() &&
  requireNamespace("rstudioapi", quietly = TRUE) &&
  rstudioapi::isAvailable()
) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

library(geomorph)
library(tidyverse)
library(janitor)
library(patchwork)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(shadowtext)

source("gm_functions.R")

#### Get Data ####
tps_files <-
  list(
    list.files(
      path = "../socal landmarking",
      pattern = "(?i)posterior.*\\.TPS$",
      full.names = TRUE
    ),
    list.files(
      path = "../matnog landmarking",
      pattern = "(?i)posterior.*\\.TPS$",
      full.names = TRUE
    )
  ) %>%
  unlist() %>%
  normalizePath()

# group label for each TPS file (same length/order as tps_files)
group <- ifelse(
  stringr::str_detect(tps_files, "socal landmarking"),
  "socal",
  "matnog"
)

# sanity checks
print(getwd())
print(table(group))
stopifnot(dir.exists("../socal landmarking"))
stopifnot(dir.exists("../matnog landmarking"))
stopifnot(length(tps_files) > 0)

#### 1. Read TPS files -> 3D array ####
Y_list <- purrr::map(
  tps_files,
  function(file) {
    Yi <- readland.tps(
      file,
      specID     = "imageID",
      negNA      = TRUE,
      readcurves = FALSE
    )
    
    # keep filename in specimen names (so we can join group labels later)
    dimnames(Yi)[[3]] <- paste0(basename(file), ":", dimnames(Yi)[[3]])
    Yi
  }
)

Y <- abind::abind(Y_list, along = 3)
print(dim(Y))  # should be 13 2 N

#### 2. Fix missing values + GPA (ONLY run GPA once) ####
Y_filled <- geomorph::estimate.missing(Y, method = "TPS")

Y.gpa <- gpagen(
  Y_filled,
  print.progress = TRUE,
  ProcD = FALSE
)

print(dim(Y.gpa$coords))  # should be 13 2 N
str(Y.gpa)

#### 3. Outline: Landmark Connections ####
outline_edges <-
  tribble(
    ~from, ~to,
    1,  5,
    5,  6,
    6,  4,
    4,  3,
    3,  2,
    2,  1,
    9,  12,
    12, 7,
    7,  8,
    8,  13,
    13, 11,
    11, 10,
    10, 9
  )

#### 3B. QUICK CHECK: mean shape socal vs matnog (overlay) ####
coords <- Y.gpa$coords

mean_socal  <- geomorph::mshape(coords[, , group == "socal",  drop = FALSE])
mean_matnog <- geomorph::mshape(coords[, , group == "matnog", drop = FALSE])

# Removed 'asp = 1' to fix the error. 
# Added color (col) so Socal is Blue and Matnog is Red.
plot(mean_socal, pch = 16, col = "blue", main = "Mean shapes: socal (blue) vs matnog (red)")
points(mean_matnog, pch = 16, col = "red")

# Add landmark numbers
text(mean_socal, labels = 1:nrow(mean_socal), pos = 3, cex = 0.8)

#### 4. Turn GPA coordinates into a tidy tibble + attach group ####
coords_tbl <- gm_step4_tidy_coords(Y.gpa$coords)

# robust mapping: specimen names (from Y.gpa) -> group (from tps_files order)
spec_to_group <- tibble(
  specimen = dimnames(Y.gpa$coords)[[3]],
  group = group
)

coords_tbl <- coords_tbl %>%
  left_join(spec_to_group, by = "specimen")

# sanity check: you should see both groups
print(table(coords_tbl$group))

#### 5. Mean shape plots ####
# (A) Your original overall mean (all specimens combined)
mean_shape_results <- gm_step5_mean_shape_plot(coords_tbl, outline_edges)

mean_shape_plot_all <- mean_shape_results$plot +
  scale_x_reverse() +
  scale_y_reverse()

mean_shape_plot_all

# (B) Compare socal vs matnog using your existing group-plot function
# If gm_step6_mean_shape_by_group() expects site_id, map group -> site_id
coords_tbl_for_groupplot <- coords_tbl %>%
  mutate(site_id = group)

mean_shape_by_group_results <- gm_step6_mean_shape_by_group(coords_tbl_for_groupplot, outline_edges)

mean_shape_plot_by_group <- mean_shape_by_group_results$plot +
  scale_x_reverse() +
  scale_y_reverse()

mean_shape_plot_by_group


#### 7. Overlay mean shapes (REVISED) ####
library(shadowtext)

coords <- Y.gpa$coords

mean_socal  <- geomorph::mshape(coords[, , group == "socal",  drop = FALSE])
mean_matnog <- geomorph::mshape(coords[, , group == "matnog", drop = FALSE])

as_lm_tbl <- function(M, grp) {
  tibble(
    landmark = seq_len(nrow(M)),
    x = M[, 1],
    y = M[, 2],
    group = grp
  )
}

mean_tbl <- bind_rows(
  as_lm_tbl(mean_matnog, "matnog"),
  as_lm_tbl(mean_socal,  "socal")
)

# --- DIRECT TRANSFORMATION ---
# To put 1 and 9 on the left:
# 1. New Y = Old X (This "stands it up")
# 2. New X = -Old Y (This moves the original vertical axis to the horizontal and flips it)
transform_upright <- function(df) {
  df %>%
    mutate(
      new_x = -y, 
      new_y = x
    ) %>%
    transmute(landmark, group, x = new_x, y = new_y)
}

mean_tbl_r <- transform_upright(mean_tbl)

# Re-calculate edges and vectors with the new coordinates
edge_df <- outline_edges %>%
  left_join(mean_tbl_r %>% rename(x_from = x, y_from = y),
            by = c("from" = "landmark")) %>%
  left_join(mean_tbl_r %>% rename(x_to = x, y_to = y),
            by = c("to" = "landmark", "group" = "group"))

vec_df <- mean_tbl_r %>%
  filter(group == "matnog") %>%
  select(landmark, x_matnog = x, y_matnog = y) %>%
  left_join(
    mean_tbl_r %>%
      filter(group == "socal") %>%
      select(landmark, x_socal = x, y_socal = y),
    by = "landmark"
  )

# Adjust label offsets for the new orientation
# Since the shape is now vertical, we want labels to sit 
# slightly outside the points on the X axis.
label_tbl <- mean_tbl_r %>%
  mutate(
    dx = if_else(x < 0, -0.02, 0.02),
    dy = if_else(y < 0, -0.02, 0.02)
  )

overlay_plot_final <-
  ggplot() +
  geom_segment(
    data = vec_df,
    aes(x = x_matnog, y = y_matnog, xend = x_socal, yend = y_socal),
    linewidth = 0.5,
    alpha = 0.45,
    arrow = arrow(length = unit(0.1, "inches"))
  ) +
  geom_segment(
    data = edge_df,
    aes(x = x_from, y = y_from, xend = x_to, yend = y_to, color = group),
    linewidth = 1
  ) +
  geom_point(
    data = mean_tbl_r,
    aes(x = x, y = y, color = group),
    size = 3
  ) +
  shadowtext::geom_shadowtext(
    data = label_tbl,
    aes(x = x + dx, y = y + dy, label = landmark, color = group),
    size = 4,
    fontface = "bold",
    bg.colour = "white",
    bg.r = 0.12,
    show.legend = FALSE
  ) +
  coord_equal() +
  labs(
    title = "Overlay of Mean Shapes",
    x = "Rotated X",
    y = "Rotated Y"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

overlay_plot_final



#### 8. Thin Plate Spline (TPS) Deformation Grid ####

# 1. Define and Transform the matrices
# We apply: new_x = -old_y and new_y = old_x
ref_raw <- mean_matnog
tar_raw <- mean_socal

# Transformation logic:
# Column 1 becomes negative Column 2
# Column 2 becomes Column 1
ref <- cbind(-ref_raw[, 2], ref_raw[, 1])
tar <- cbind(-tar_raw[, 2], tar_raw[, 1])

# 2. Convert edges to matrix (if not already done)
outline_edges_matrix <- as.matrix(outline_edges[, c("from", "to")])

# 3. Define grid parameters
gp <- gridPar(grid.col = "grey70", grid.lwd = 0.5, pt.bg = "black", pt.size = 1.5)

# 4. Plot the TPS grid
plotRefToTarget(ref, tar, 
                method = "TPS", 
                gridPars = gp,
                mag = 1.0, 
                links = outline_edges_matrix)

title(main = "Shape Deformation: Matnog to Socal")



#### 9. PCA Plot (High-Contrast Version) ####

# 1. Run the PCA on the GPA coordinates
pca_res <- geomorph::gm.prcomp(Y.gpa$coords)

# 2. Extract scores into a dataframe
pca_df <- tibble(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  group = group
)

# 3. Calculate % variance
pc1_var <- round(pca_res$pc.summary$importance[2, 1] * 100, 1)
pc2_var <- round(pca_res$pc.summary$importance[2, 2] * 100, 1)

# 4. Create the plot
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group, fill = group)) +
  # Background grid lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey80") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
  
  # IMPROVED ELLIPSES: Higher alpha and added a border line
  stat_ellipse(geom = "polygon", 
               alpha = 0.25,        # Made the fill darker
               linewidth = 0.8,     # Added a visible border line
               aes(fill = group, color = group), 
               linetype = "solid") + 
  
  # IMPROVED POINTS: Added a black outline for "pop"
  geom_point(size = 3.5, shape = 21, color = "black", stroke = 0.5) +
  
  # SATURATED COLORS: Stronger Red and Blue
  scale_color_manual(values = c("matnog" = "#E41A1C", "socal" = "#377EB8")) +
  scale_fill_manual(values = c("matnog" = "#E41A1C", "socal" = "#377EB8")) +
  
  labs(
    title = "Principal Component Analysis of Shape",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    panel.grid.minor = element_blank() # Cleaned up the grid
  )

pca_plot



#### 6. Mean shape by era x site with ggplot ####

mean_shape_era_site_results <- gm_step6_mean_shape_by_group(coords_tbl, outline_edges)
mean_shape_tbl_era_site <- mean_shape_era_site_results$mean_shape_tbl
outline_mean_tbl_era_site <- mean_shape_era_site_results$outline_mean_tbl
mean_shape_plot_era_site <- mean_shape_era_site_results$plot



#### 7. PCA with Manual Grouping (FIXED) ####

# 1. Run the PCA math
pca_results <- gm_step7_pca(Y.gpa$coords, gm_parse_specimen)
pca <- pca_results$pca
scores_tbl <- pca_results$scores_tbl

# 2. OVERWRITE the site_id based on your folder structure
# We look at the tps_files list we created in the "Get Data" step
scores_tbl <- scores_tbl %>%
  mutate(site_id = if_else(grepl("matnog", tps_files), "matnog", "socal"))

# 3. Plot with the corrected site_id
pca_plot_final <- ggplot(scores_tbl, aes(x = PC1, y = PC2, color = site_id, fill = site_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey80") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
  
  # Now this will see 'matnog' and 'socal' and draw TWO circles
  stat_ellipse(geom = "polygon", alpha = 0.25, linewidth = 0.8) + 
  
  geom_point(size = 3.5, shape = 21, color = "black", stroke = 0.5) +
  
  # Apply your site colors
  scale_color_manual(values = c("matnog" = "#E41A1C", "socal" = "#377EB8")) +
  scale_fill_manual(values = c("matnog" = "#E41A1C", "socal" = "#377EB8")) +
  
  labs(
    title = "PCA of Vertebra Shape",
    subtitle = "Comparing Matnog vs. Socal Sites",
    x = "PC1", y = "PC2"
  ) +
  theme_minimal(base_size = 14)

pca_plot_final
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
