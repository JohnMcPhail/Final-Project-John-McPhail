# Utility helpers and step-wise functions for vertebra landmark workflows.
# Steps correspond to process_landmarks_vertebra_2nd-to-last.R sections.
# Requires tidyverse, geomorph, and patchwork to be attached by the caller.

gm_parse_specimen <-
  function(data)
  {
    data %>%
      mutate(
        file_name =
          str_replace(
            specimen,
            "\\.TPS:sde_xray.*",
            "\\.TPS"
          ),
        individual_id =
          str_remove(
            specimen,
            "^.*\\.TPS:"
          ) %>%
          str_remove(
            .,
            "^sde_xray_...\\-._Sde\\-"
          ) %>%
          str_remove(
            .,
            "^sde_xray_...\\..\\-._Sde\\-"
          ) %>%
          str_remove(
            .,
            "\\.png"
          ),
        era = str_sub(
          individual_id,
          1,
          1
        ),
        site_id = str_sub(
          individual_id,
          2,
          4
        )
      )
  }

gm_outline_edges <-
  function()
  {
    tribble(
      ~from, ~to,
      7, 5,
      5, 1,
      1, 2,
      2, 3,
      3, 4,
      4, 6,
      6, 8,
      15, 10,
      10, 9,
      9, 13,
      13, 14,
      14, 12,
      12, 11,
      11, 16
    )
  }

# Step 4: Turn GPA coordinates into a tidy tibble
gm_step4_tidy_coords <-
  function(
    coords_array,
    parse_specimen_fn = gm_parse_specimen
  )
  {
    coords_array %>%
      as.data.frame.table(responseName = "value") %>%
      as_tibble() %>%
      mutate(
        landmark = as.integer(Var1),
        axis     = if_else(Var2 == "X", "y", "x"),
        specimen = as.character(Var3)
      ) %>%
      parse_specimen_fn() %>%
      select(
        landmark,
        axis,
        file_name,
        era,
        site_id,
        individual_id,
        specimen,
        value
      ) %>%
      pivot_wider(
        names_from  = axis,
        values_from = value
      ) %>%
      mutate(y = -y)
  }

# Step 5: Mean shape + all specimens plot
gm_step5_mean_shape_plot <-
  function(
    coords_tbl,
    outline_edges = gm_outline_edges(),
    title = "Mean shape (red) with all aligned specimens"
  )
  {
    mean_shape_tbl <-
      coords_tbl %>%
      group_by(landmark) %>%
      summarise(
        x = mean(x, na.rm = TRUE),
        y = mean(y, na.rm = TRUE),
        .groups = "drop"
      )
    
    outline_mean_tbl <-
      outline_edges %>%
      left_join(
        mean_shape_tbl,
        by = c("from" = "landmark")
      ) %>%
      rename(
        x_from = x,
        y_from = y
      ) %>%
      left_join(
        mean_shape_tbl,
        by = c("to" = "landmark")
      ) %>%
      rename(
        x_to = x,
        y_to = y
      )
    
    plot_mean_shape <- ggplot() +
      geom_point(
        data = coords_tbl,
        aes(x = x, y = y),
        alpha = 0.15,
        size = 1
      ) +
      geom_segment(
        data = outline_mean_tbl,
        aes(
          x = x_from,
          y = y_from,
          xend = x_to,
          yend = y_to
        ),
        linewidth = 0.6,
        color = "red"
      ) +
      geom_point(
        data = mean_shape_tbl,
        aes(x = x, y = y),
        color = "red",
        size = 2
      ) +
      geom_text(
        data = mean_shape_tbl,
        aes(
          x = x,
          y = y,
          label = landmark
        ),
        vjust = -0.5,
        size = 3
      ) +
      coord_equal() +
      labs(
        x = "X (Procrustes)",
        y = "Y (Procrustes)",
        title = title
      ) +
      theme_minimal() +
      facet_grid(era ~ site_id)
    
    list(
      mean_shape_tbl = mean_shape_tbl,
      outline_mean_tbl = outline_mean_tbl,
      plot = plot_mean_shape
    )
  }

# Step 6: Mean shape by era x site plot
gm_step6_mean_shape_by_group <-
  function(
    coords_tbl,
    outline_edges = gm_outline_edges()
  )
  {
    mean_shape_tbl_era_site <-
      coords_tbl %>%
      group_by(era, site_id, landmark) %>%
      summarise(
        x = mean(x, na.rm = TRUE),
        y = mean(y, na.rm = TRUE),
        .groups = "drop"
      )
    
    outline_mean_tbl_era_site <-
      outline_edges %>%
      tidyr::expand_grid(
        era     = unique(coords_tbl$era),
        site_id = unique(coords_tbl$site_id)
      ) %>%
      left_join(
        mean_shape_tbl_era_site,
        by = c(
          "era",
          "site_id",
          "from" = "landmark"
        )
      ) %>%
      rename(
        x_from = x,
        y_from = y
      ) %>%
      left_join(
        mean_shape_tbl_era_site,
        by = c(
          "era",
          "site_id",
          "to" = "landmark"
        )
      ) %>%
      rename(
        x_to = x,
        y_to = y
      )
    
    plot_mean_shape_era_site <- ggplot() +
      geom_point(
        data = coords_tbl,
        aes(
          x = x,
          y = y,
          color = site_id
        ),
        alpha = 0.15,
        size = 1
      ) +
      geom_segment(
        data = outline_mean_tbl_era_site,
        aes(
          x = x_from,
          y = y_from,
          xend = x_to,
          yend = y_to,
          color = site_id,
          group = interaction(era, site_id)
        ),
        linewidth = 0.6
      ) +
      geom_point(
        data = mean_shape_tbl_era_site,
        aes(
          x = x,
          y = y,
          color = site_id
        ),
        color = "red",
        size = 2
      ) +
      geom_text(
        data = mean_shape_tbl_era_site,
        aes(
          x = x,
          y = y,
          label = landmark,
          color = site_id
        ),
        vjust = -0.5,
        size = 3
      ) +
      coord_equal() +
      labs(
        x = "X (Procrustes)",
        y = "Y (Procrustes)",
        title = "Mean shapes by era and site"
      ) +
      theme_minimal() +
      facet_grid(era ~ site_id)
    
    list(
      mean_shape_tbl = mean_shape_tbl_era_site,
      outline_mean_tbl = outline_mean_tbl_era_site,
      plot = plot_mean_shape_era_site
    )
  }

# Step 7: PCA with gm.prcomp + tidyverse
gm_step7_pca <-
  function(
    coords_array,
    parse_specimen_fn = gm_parse_specimen,
    ellipse_level = 0.68,
    title = "PCA of vertebra shape (PC1 vs PC2)"
  )
  {
    pca <- gm.prcomp(coords_array)
    
    scores_tbl <-
      pca$x %>%
      as.data.frame() %>%
      rownames_to_column("specimen") %>%
      as_tibble() %>%
      parse_specimen_fn() %>%
      rename_with(~ str_replace(.x, "^Comp", "PC")) %>%
      mutate(era_site = interaction(era, site_id, drop = TRUE))
    
    plot_pca <-
      scores_tbl %>%
      ggplot() +
      aes(
        x = PC1,
        y = PC2,
        color = site_id,
        shape = era
      ) +
      stat_ellipse(
        aes(group = era_site),
        type = "norm",
        level = ellipse_level,
        linewidth = 0.7
      ) +
      geom_point(
        size = 2,
        alpha = 0.8
      ) +
      coord_equal() +
      labs(
        x = "PC1",
        y = "PC2",
        title = title
      ) +
      theme_minimal()
    
    plot_pca_34 <-
      scores_tbl %>%
      ggplot() +
      aes(
        x = PC3,
        y = PC4,
        color = site_id,
        shape = era
      ) +
      geom_point(
        size = 2,
        alpha = 0.8
      ) +
      coord_equal() +
      labs(
        x = "PC3",
        y = "PC4",
        title = "PCA of vertebra shape (PC3 vs PC4)"
      ) +
      theme_minimal()
    
    list(
      pca = pca,
      scores_tbl = scores_tbl,
      plot_pca = plot_pca,
      plot_pca_34 = plot_pca_34
    )
  }

gm_make_shape_tbl <-
  function(shape_mat)
  {
    as.data.frame(shape_mat) %>%
      as_tibble(.name_repair = "minimal") %>%
      mutate(
        landmark = row_number(),
        x = Y,
        y = -X
      ) %>%
      select(landmark, x, y)
  }

gm_make_outline_tbl <-
  function(
    shape_tbl,
    outline_edges
  )
  {
    outline_edges %>%
      left_join(
        shape_tbl,
        by = c("from" = "landmark")
      ) %>%
      rename(
        x_from = x,
        y_from = y
      ) %>%
      left_join(
        shape_tbl,
        by = c("to" = "landmark")
      ) %>%
      rename(
        x_to = x,
        y_to = y
      )
  }

gm_shape_plot <-
  function(
    shape_mat,
    title = "",
    outline_edges,
    mean_shape_tbl,
    outline_mean_tbl
  )
  {
    shape_tbl <- gm_make_shape_tbl(shape_mat)
    outline_tbl <- gm_make_outline_tbl(shape_tbl, outline_edges)
    
    ggplot() +
      geom_segment(
        data = outline_mean_tbl,
        aes(
          x = x_from,
          y = y_from,
          xend = x_to,
          yend = y_to
        ),
        linewidth = 0.6,
        color = "grey70"
      ) +
      geom_point(
        data = mean_shape_tbl,
        aes(x = x, y = y),
        color = "grey60",
        size = 1.5
      ) +
      geom_segment(
        data = outline_tbl,
        aes(
          x = x_from,
          y = y_from,
          xend = x_to,
          yend = y_to
        )
      ) +
      geom_point(
        data = shape_tbl,
        aes(x = x, y = y),
        size = 2
      ) +
      coord_equal() +
      ggtitle(title) +
      theme_void()
  }

# Step 8: Make shapes at ends of PC axes
gm_step8_pc_end_shapes <-
  function(
    pca,
    outline_edges,
    mean_shape_tbl,
    outline_mean_tbl
  )
  {
    list(
      shape_pc1_min = gm_shape_plot(
        pca$shapes$shapes.comp1$min,
        "PC1 low",
        outline_edges,
        mean_shape_tbl,
        outline_mean_tbl
      ),
      shape_pc1_max = gm_shape_plot(
        pca$shapes$shapes.comp1$max,
        "PC1 high",
        outline_edges,
        mean_shape_tbl,
        outline_mean_tbl
      ),
      shape_pc2_min = gm_shape_plot(
        pca$shapes$shapes.comp2$min,
        "PC2 low",
        outline_edges,
        mean_shape_tbl,
        outline_mean_tbl
      ),
      shape_pc2_max = gm_shape_plot(
        pca$shapes$shapes.comp2$max,
        "PC2 high",
        outline_edges,
        mean_shape_tbl,
        outline_mean_tbl
      )
    )
  }

# Step 9: Corner shapes in PC1–PC2 space
gm_step9_corner_shapes <-
  function(
    pca,
    consensus_mat,
    outline_edges,
    mean_shape_tbl,
    outline_mean_tbl
  )
  {
    d1_min <- pca$shapes$shapes.comp1$min - consensus_mat
    d1_max <- pca$shapes$shapes.comp1$max - consensus_mat
    d2_min <- pca$shapes$shapes.comp2$min - consensus_mat
    d2_max <- pca$shapes$shapes.comp2$max - consensus_mat
    
    shape_pc1min_pc2min <- consensus_mat + d1_min + d2_min
    shape_pc1min_pc2max <- consensus_mat + d1_min + d2_max
    shape_pc1max_pc2min <- consensus_mat + d1_max + d2_min
    shape_pc1max_pc2max <- consensus_mat + d1_max + d2_max
    
    list(
      shape_corner_ll = gm_shape_plot(
        shape_pc1min_pc2min,
        "PC1 low, PC2 low",
        outline_edges,
        mean_shape_tbl,
        outline_mean_tbl
      ),
      shape_corner_lu = gm_shape_plot(
        shape_pc1min_pc2max,
        "PC1 low, PC2 high",
        outline_edges,
        mean_shape_tbl,
        outline_mean_tbl
      ),
      shape_corner_rl = gm_shape_plot(
        shape_pc1max_pc2min,
        "PC1 high, PC2 low",
        outline_edges,
        mean_shape_tbl,
        outline_mean_tbl
      ),
      shape_corner_ru = gm_shape_plot(
        shape_pc1max_pc2max,
        "PC1 high, PC2 high",
        outline_edges,
        mean_shape_tbl,
        outline_mean_tbl
      )
    )
  }

# Step 10: Arrange shapes around the PCA plot
gm_step10_arrange_pca_shapes <-
  function(
    plot_pca,
    shapes_pc,
    corner_shapes
  )
  {
    axis_layout <- (
      plot_spacer() | shapes_pc$shape_pc2_max | plot_spacer()
    ) /
      (
        shapes_pc$shape_pc1_min | plot_pca | shapes_pc$shape_pc1_max
      ) /
      (
        plot_spacer() | shapes_pc$shape_pc2_min | plot_spacer()
      )
    
    corner_layout <- (
      corner_shapes$shape_corner_lu | plot_spacer() | corner_shapes$shape_corner_ru
    ) /
      (
        plot_spacer() | plot_pca | plot_spacer()
      ) /
      (
        corner_shapes$shape_corner_ll | plot_spacer() | corner_shapes$shape_corner_rl
      ) +
      plot_layout(
        widths  = c(1, 1, 1),
        heights = c(1, 1, 1)
      )
    
    list(
      axis_layout = axis_layout,
      corner_layout = corner_layout
    )
  }

# Step 11.1: Build a geomorph data frame
gm_step11_1_geomorph_df <-
  function(
    coords_array,
    Csize,
    coords_tbl
  )
  {
    specimen_order <- dimnames(coords_array)[[3]]
    if (is.null(specimen_order)) {
      stop("coords_array needs specimen names in dimnames(coords_array)[[3]]")
    }
    
    spec_meta <-
      coords_tbl %>%
      distinct(specimen, era, site_id) %>%
      mutate(specimen = factor(specimen, levels = specimen_order)) %>%
      arrange(specimen)
    
    stopifnot(all(as.character(spec_meta$specimen) == specimen_order))
    
    gdf <-
      geomorph.data.frame(
        coords  = coords_array,
        Csize   = Csize,
        era     = spec_meta$era,
        site_id = spec_meta$site_id
      )
    
    list(
      gdf = gdf,
      specimen_order = specimen_order,
      spec_meta = spec_meta
    )
  }

# Step 11.2: Test for overall allometry
gm_step11_2_test_allometry <-
  function(
    gdf,
    iter = 9999
  )
  {
    procD.lm(
      coords ~ log(Csize),
      data  = gdf,
      iter  = iter
    )
  }

# Step 11.3: Does allometry differ by era or site?
gm_step11_3_test_allometry_by_group <- 
  function(
    gdf,
    iter = 9999
  ) {
    fit_allom_era <-
      procD.lm(
        coords ~ log(Csize) * era,
        data  = gdf,
        iter  = iter
      )
    
    fit_allom_site <-
      procD.lm(
        coords ~ log(Csize) * site_id,
        data  = gdf,
        iter  = iter
      )
    
    fit_allom_full <-
      procD.lm(
        coords ~ log(Csize) * era * site_id,
        data  = gdf,
        iter  = iter
      )
    
    list(
      fit_allom_era = fit_allom_era,
      fit_allom_site = fit_allom_site,
      fit_allom_full = fit_allom_full
    )
  }

# Step 11.4: Visualizing the allometric pattern
gm_step11_4_allometry_plot <- 
  function(
    fit_allom,
    gdf
  ) {
    plot_obj <- plotAllometry(
      fit_allom,
      size   = gdf$Csize,
      logsz  = TRUE,
      method = "RegScore"
    )
    
    allom_df <- tibble(
      logCsize    = log(gdf$Csize),
      RegScore    = as.numeric(plot_obj$RegScore),
      era         = gdf$era,
      site_id     = gdf$site_id
    )
    
    allom_plot <-
      ggplot(
        allom_df,
        aes(
          x = logCsize,
          y = RegScore,
          color = era,
          shape = site_id
        )
      ) +
      geom_point(size = 2) +
      geom_smooth(method = "lm", se = FALSE) +
      labs(
        x = "log(Centroid Size)",
        y = "Shape regression score",
        color = "Era",
        title = "Allometry of vertebra shape by era"
      ) +
      theme_minimal()
    
    list(
      reg_scores = plot_obj$RegScore,
      allom_df = allom_df,
      plot = allom_plot
    )
  }

# Step 11.5: size-corrected shape dataset
gm_step11_5_size_corrected_shapes <- 
  function(
    fit_allom,
    coords_array,
    reference_logCsize = NULL
  ) {
    res_mat <- residuals(fit_allom)
    
    p <- dim(coords_array)[1]
    k <- dim(coords_array)[2]
    
    if (is.null(reference_logCsize)) {
      if (!is.null(fit_allom$data$Csize)) {
        reference_logCsize <-
          fit_allom$data$Csize %>%
          log() %>%
          mean(na.rm = TRUE)
      } else if (
        is.matrix(fit_allom$X) &&
          "log(Csize)" %in% colnames(fit_allom$X)
      ) {
        reference_logCsize <-
          fit_allom$X[, "log(Csize)"] %>%
          mean(na.rm = TRUE)
      } else {
        reference_logCsize <- NA_real_
      }
    }
    
    coef_mat <- fit_allom$coefficients
    
    ref_shape_vec <-
      if (
        is.matrix(coef_mat) &&
          is.finite(reference_logCsize) &&
          "(Intercept)" %in% rownames(coef_mat) &&
          "log(Csize)" %in% rownames(coef_mat)
      ) {
        coef_mat["(Intercept)", ] +
          coef_mat["log(Csize)", ] * reference_logCsize
      } else {
        if (is.null(fit_allom$fitted)) {
          stop("fit_allom$fitted is missing; cannot build size-adjusted shapes")
        }
        
        colMeans(fit_allom$fitted)
      }
    
    size_adjusted_mat <-
      res_mat +
      matrix(
        ref_shape_vec,
        nrow = nrow(res_mat),
        ncol = ncol(res_mat),
        byrow = TRUE
      )
    
    size_adjusted_shapes <- arrayspecs(size_adjusted_mat, p = p, k = k)
    dimnames(size_adjusted_shapes) <- dimnames(coords_array)
    
    attr(
      size_adjusted_shapes,
      "reference_logCsize"
    ) <- reference_logCsize
    
    size_adjusted_shapes
  }

# Step 11.6: Test era and site on size-corrected shapes
gm_step11_6_test_residual_shapes <- 
  function(
    shape_resid,
    spec_meta,
    iter = 9999
  ) {
    specimen_order <- dimnames(shape_resid)[[3]]
    if (is.null(specimen_order)) {
      stop("shape_resid needs specimen names in dimnames(shape_resid)[[3]]")
    }
    
    spec_meta <- spec_meta %>%
      mutate(specimen = factor(specimen, levels = specimen_order)) %>%
      arrange(specimen)
    
    stopifnot(all(as.character(spec_meta$specimen) == specimen_order))
    
    gdf_resid <- geomorph.data.frame(
      coords  = shape_resid,
      era     = spec_meta$era,
      site_id = spec_meta$site_id
    )
    
    fit_resid_era_site <- procD.lm(
      coords ~ era * site_id,
      data  = gdf_resid,
      iter  = iter
    )
    
    list(
      gdf_resid = gdf_resid,
      fit_resid_era_site = fit_resid_era_site
    )
  }

