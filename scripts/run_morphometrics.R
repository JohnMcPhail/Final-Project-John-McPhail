#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("geomorph", quietly = TRUE)) {
    stop(
      "Missing required package: geomorph\n",
      "Install in R with: install.packages('geomorph')\n",
      call. = FALSE
    )
  }
  library(geomorph)
})

default_matnog_dir <- file.path(
  "data",
  "raw",
  "johns_landmarked_photos",
  "2024-12_nmnh_c-matnog_a-matnog-sorsogon_xrays",
  "xrays_cropped-for-geomorph"
)

default_sacol_dir <- file.path(
  "data",
  "raw",
  "johns_landmarked_photos",
  "2024-12_nmnh_c-taluksangay_a-sacol-island_xrays",
  "xrays_cropped-for-geomorph"
)

parse_args <- function(args) {
  opts <- list(
    matnog_dir = default_matnog_dir,
    sacol_dir = default_sacol_dir,
    out_dir = "results",
    processed_dir = file.path("data", "processed"),
    pattern = "Back\\.TPS$",
    recursive = FALSE,
    iter = 999,
    seed = 1L
  )

  for (arg in args) {
    if (arg %in% c("-h", "--help")) {
      cat(
        "Usage: Rscript scripts/run_morphometrics.R [options]\n\n",
        "Options:\n",
        "  --matnog-dir=PATH      Input dir for Matnog TPS files\n",
        "  --sacol-dir=PATH       Input dir for Sacol TPS files\n",
        "  --pattern=REGEX        File match pattern (default: Back\\\\.TPS$)\n",
        "  --recursive=TRUE|FALSE Search for files recursively (default: FALSE)\n",
        "  --out-dir=PATH         Results directory (default: results)\n",
        "  --processed-dir=PATH   Processed data directory (default: data/processed)\n",
        "  --iter=N               Permutations for procD.lm (default: 999)\n",
        "  --seed=N               RNG seed (default: 1)\n",
        sep = ""
      )
      quit(status = 0)
    }

    if (startsWith(arg, "--matnog-dir=")) opts$matnog_dir <- sub("^--matnog-dir=", "", arg)
    if (startsWith(arg, "--sacol-dir=")) opts$sacol_dir <- sub("^--sacol-dir=", "", arg)
    if (startsWith(arg, "--pattern=")) opts$pattern <- sub("^--pattern=", "", arg)
    if (startsWith(arg, "--out-dir=")) opts$out_dir <- sub("^--out-dir=", "", arg)
    if (startsWith(arg, "--processed-dir=")) opts$processed_dir <- sub("^--processed-dir=", "", arg)

    if (startsWith(arg, "--recursive=")) {
      val <- tolower(sub("^--recursive=", "", arg))
      opts$recursive <- isTRUE(val %in% c("true", "t", "1", "yes", "y"))
    }

    if (startsWith(arg, "--iter=")) opts$iter <- as.integer(sub("^--iter=", "", arg))
    if (startsWith(arg, "--seed=")) opts$seed <- as.integer(sub("^--seed=", "", arg))
  }

  opts
}

get_tps_files <- function(dir, pattern, recursive) {
  if (!dir.exists(dir)) {
    stop("Input directory not found: ", dir, call. = FALSE)
  }
  list.files(dir, pattern = pattern, full.names = TRUE, recursive = recursive)
}

read_tps_landmarks <- function(files) {
  tps_list <- lapply(
    files,
    function(path) {
      suppressWarnings(utils::capture.output(
        res <- geomorph::readland.tps(path, specID = "None", warnmsg = FALSE)
      ))
      res
    }
  )
  dims <- vapply(tps_list, function(x) paste(dim(x), collapse = "x"), character(1))
  if (length(unique(dims)) != 1) {
    stop(
      "TPS files contain inconsistent landmark dimensions.\n",
      "Found: ", paste(unique(dims), collapse = ", "),
      call. = FALSE
    )
  }

  n_spec <- length(tps_list)
  n_land <- nrow(tps_list[[1]])
  k_dim <- ncol(tps_list[[1]])

  all_landmarks <- array(NA_real_, dim = c(n_land, k_dim, n_spec))
  for (i in seq_len(n_spec)) all_landmarks[, , i] <- tps_list[[i]]

  list(all_landmarks = all_landmarks, n_land = n_land, k_dim = k_dim)
}

save_png <- function(path, plot_fun, width = 2000, height = 1600, res = 200) {
  grDevices::png(filename = path, width = width, height = height, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot_fun()
}

draw_ellipse <- function(center, covmat, level = 0.95, npoints = 200, ...) {
  theta <- seq(0, 2 * pi, length.out = npoints)
  circle <- cbind(cos(theta), sin(theta))
  eig <- eigen(covmat)
  radii <- sqrt(eig$values * stats::qchisq(level, df = 2))
  shape <- circle %*% diag(radii) %*% t(eig$vectors)
  graphics::lines(center[1] + shape[, 1], center[2] + shape[, 2], ...)
}

get_anova_value <- function(anova_tbl, term, col) {
  tbl <- anova_tbl
  if (is.list(tbl) && !is.null(tbl$table)) tbl <- tbl$table
  if (is.null(rownames(tbl)) || is.null(colnames(tbl))) return(NA_real_)
  if (!term %in% rownames(tbl)) return(NA_real_)
  if (!col %in% colnames(tbl)) return(NA_real_)
  as.numeric(tbl[term, col])
}

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.4f", p)
}

run_analysis <- function(opts) {
  set.seed(opts$seed)

  matnog_files <- get_tps_files(opts$matnog_dir, opts$pattern, opts$recursive)
  sacol_files <- get_tps_files(opts$sacol_dir, opts$pattern, opts$recursive)

  if (length(matnog_files) == 0) stop("No TPS files found in: ", opts$matnog_dir, call. = FALSE)
  if (length(sacol_files) == 0) stop("No TPS files found in: ", opts$sacol_dir, call. = FALSE)

  tps_files <- c(matnog_files, sacol_files)
  group <- factor(
    c(rep("Matnog", length(matnog_files)), rep("Sacol", length(sacol_files))),
    levels = c("Matnog", "Sacol")
  )

  read_res <- read_tps_landmarks(tps_files)
  all_landmarks <- read_res$all_landmarks
  n_land <- read_res$n_land
  k_dim <- read_res$k_dim

  gpa <- geomorph::gpagen(all_landmarks, print.progress = FALSE)
  coords <- gpa$coords
  CS <- gpa$Csize

  model_data <- list(coords = coords, CS = CS, group = group)

  pca <- geomorph::gm.prcomp(coords)
  scores <- pca$x
  var_expl <- 100 * (pca$d ^ 2 / sum(pca$d ^ 2))

  mean_matnog <- apply(coords[, , group == "Matnog", drop = FALSE], c(1, 2), mean)
  mean_sacol <- apply(coords[, , group == "Sacol", drop = FALSE], c(1, 2), mean)
  line_order <- seq_len(nrow(mean_matnog))

  procdist_to_mean <- numeric(dim(coords)[3])
  for (i in seq_len(dim(coords)[3])) {
    mean_shape_i <- if (group[i] == "Matnog") mean_matnog else mean_sacol
    procdist_to_mean[i] <- sqrt(sum((coords[, , i] - mean_shape_i) ^ 2))
  }

  allom_fit <- geomorph::procD.lm(coords ~ CS, data = model_data, iter = opts$iter, print.progress = FALSE)
  shape_resid_array <- geomorph::arrayspecs(allom_fit$residuals, p = n_land, k = k_dim)

  mean_matnog_resid <- apply(shape_resid_array[, , group == "Matnog", drop = FALSE], c(1, 2), mean)
  mean_sacol_resid <- apply(shape_resid_array[, , group == "Sacol", drop = FALSE], c(1, 2), mean)

  coords_means <- array(NA_real_, dim = c(n_land, k_dim, 2))
  coords_means[, , 1] <- mean_matnog_resid
  coords_means[, , 2] <- mean_sacol_resid
  gpa_means <- geomorph::gpagen(coords_means, print.progress = FALSE)
  mean_matnog_resid_al <- gpa_means$coords[, , 1]
  mean_sacol_resid_al <- gpa_means$coords[, , 2]

  dir.create(opts$out_dir, showWarnings = FALSE, recursive = TRUE)
  fig_dir <- file.path(opts$out_dir, "figures")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(opts$processed_dir, showWarnings = FALSE, recursive = TRUE)

  specimen_id <- sub("_Back$", "", tools::file_path_sans_ext(basename(tps_files)))
  specimen_tbl <- data.frame(
    specimen_id = specimen_id,
    group = group,
    tps_path = tps_files,
    stringsAsFactors = FALSE
  )

  utils::write.csv(specimen_tbl, file.path(opts$processed_dir, "specimens.csv"), row.names = FALSE)
  saveRDS(
    list(
      tps_files = tps_files,
      group = group,
      coords_raw = all_landmarks,
      coords_gpa = coords,
      centroid_size = CS,
      pca_scores = scores
    ),
    file.path(opts$processed_dir, "geomorph_data.rds")
  )

  save_png(file.path(fig_dir, "pca_pc1_pc2.png"), plot_fun = function() {
    cols <- c("Matnog" = "blue", "Sacol" = "red")
    pch_vals <- c("Matnog" = 16, "Sacol" = 17)

    pad <- 0.50
    xrange <- range(scores[, 1])
    yrange <- range(scores[, 2])
    xlim_zoom <- xrange + c(-pad, pad) * diff(xrange)
    ylim_zoom <- yrange + c(-pad, pad) * diff(yrange)

    graphics::plot(
      scores[, 1],
      scores[, 2],
      col = cols[group],
      pch = pch_vals[group],
      xlim = xlim_zoom,
      ylim = ylim_zoom,
      xlab = paste0("PC1 (", round(var_expl[1], 1), "%)"),
      ylab = paste0("PC2 (", round(var_expl[2], 1), "%)"),
      main = "PCA of vertebra landmark shape\nMatnog vs Sacol"
    )

    graphics::legend(
      "topright",
      legend = levels(group),
      col = cols[levels(group)],
      pch = pch_vals[levels(group)],
      bty = "n"
    )

    for (g in levels(group)) {
      idx <- which(group == g)
      center <- colMeans(scores[idx, 1:2, drop = FALSE])
      covmat <- stats::cov(scores[idx, 1:2, drop = FALSE])
      draw_ellipse(center, covmat, col = cols[g], lwd = 2)
    }

    graphics::text(scores[, 1], scores[, 2], labels = specimen_id, pos = 3, cex = 0.55)
  })

  save_png(file.path(fig_dir, "mean_shape_raw.png"), plot_fun = function() {
    graphics::plot(
      mean_matnog,
      asp = 1,
      type = "n",
      xlab = "X",
      ylab = "Y",
      main = "Mean vertebra shape (raw)\nMatnog vs Sacol"
    )
    graphics::lines(mean_matnog[c(line_order, line_order[1]), ], col = "blue", lwd = 2)
    graphics::points(mean_matnog[line_order, ], col = "blue", pch = 16)
    graphics::lines(mean_sacol[c(line_order, line_order[1]), ], col = "red", lwd = 2, lty = 2)
    graphics::points(mean_sacol[line_order, ], col = "red", pch = 17)
    graphics::legend(
      "topright",
      legend = c("Matnog mean", "Sacol mean"),
      col = c("blue", "red"),
      lty = c(1, 2),
      pch = c(16, 17),
      bty = "n"
    )
  })

  save_png(file.path(fig_dir, "within_group_variation.png"), plot_fun = function() {
    spread_df <- data.frame(group = group, dist = procdist_to_mean)
    graphics::boxplot(
      dist ~ group,
      data = spread_df,
      xlab = "Population",
      ylab = "Procrustes distance to group mean",
      main = "Within-group shape variation"
    )
  })

  save_png(file.path(fig_dir, "mean_shape_size_corrected.png"), plot_fun = function() {
    line_order2 <- seq_len(nrow(mean_matnog_resid_al))
    pad <- 0.015
    all_coords <- rbind(mean_matnog_resid_al, mean_sacol_resid_al)
    xrange <- range(all_coords[, 1])
    yrange <- range(all_coords[, 2])
    xlim_zoom <- xrange + c(-pad, pad) * diff(xrange)
    ylim_zoom <- yrange + c(-pad, pad) * diff(yrange)

    graphics::plot(
      mean_matnog_resid_al,
      asp = 1,
      type = "n",
      xlim = xlim_zoom,
      ylim = ylim_zoom,
      xlab = "X coordinate (allometry corrected)",
      ylab = "Y coordinate (allometry corrected)",
      main = "Size-corrected mean vertebra shape\nMatnog vs Sacol"
    )

    graphics::lines(mean_matnog_resid_al[c(line_order2, line_order2[1]), ], col = "blue", lwd = 2)
    graphics::points(mean_matnog_resid_al[line_order2, ], col = "blue", pch = 16)

    graphics::lines(mean_sacol_resid_al[c(line_order2, line_order2[1]), ], col = "red", lwd = 2, lty = 2)
    graphics::points(mean_sacol_resid_al[line_order2, ], col = "red", pch = 17)

    graphics::legend(
      "topright",
      legend = c("Matnog (size-corrected)", "Sacol (size-corrected)"),
      col = c("blue", "red"),
      lty = c(1, 2),
      pch = c(16, 17),
      bty = "n"
    )
  })

  shape_group_fit <- geomorph::procD.lm(coords ~ group, data = model_data, iter = opts$iter, print.progress = FALSE)
  shape_group_anova <- anova(shape_group_fit)

  size_fit <- stats::lm(CS ~ group)
  size_r2 <- summary(size_fit)$r.squared
  size_p <- stats::anova(size_fit)[["Pr(>F)"]][1]

  shape_group_size_fit <- geomorph::procD.lm(coords ~ group + CS, data = model_data, iter = opts$iter, print.progress = FALSE)
  shape_group_size_anova <- anova(shape_group_size_fit)

  summary_lines <- c(
    "================= SUMMARY =================",
    "",
    sprintf("TPS files: Matnog=%d, Sacol=%d (total=%d)", length(matnog_files), length(sacol_files), length(tps_files)),
    "",
    "1) Shape difference between Matnog and Sacol",
    sprintf("   - Model: coords ~ group (perm=%d)", opts$iter),
    sprintf("   - R² = %0.4f", get_anova_value(shape_group_anova, "group", "Rsq")),
    sprintf("   - p = %s", fmt_p(get_anova_value(shape_group_anova, "group", "Pr(>F)"))),
    "",
    "2) Centroid size difference between groups",
    "   - Model: CS ~ group",
    sprintf("   - R² = %0.4f", size_r2),
    sprintf("   - p = %s", fmt_p(size_p)),
    "",
    "3) Shape differences after correcting for size (allometry)",
    sprintf("   - Model: coords ~ group + CS (perm=%d)", opts$iter),
    sprintf("   - group p = %s", fmt_p(get_anova_value(shape_group_size_anova, "group", "Pr(>F)"))),
    sprintf("   - CS    p = %s", fmt_p(get_anova_value(shape_group_size_anova, "CS", "Pr(>F)"))),
    "",
    "Outputs:",
    sprintf("  - Figures: %s", fig_dir),
    sprintf("  - Processed: %s", opts$processed_dir),
    "==========================================="
  )

  writeLines(summary_lines, con = file.path(opts$out_dir, "summary.txt"))
  writeLines(summary_lines)

  invisible(list(
    coords = coords,
    group = group,
    centroid_size = CS,
    pca = pca,
    anova_shape_group = shape_group_anova,
    anova_shape_group_size = shape_group_size_anova,
    size_model = size_fit
  ))
}

is_executed_as_script <- function() {
  file_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_args) == 0) return(FALSE)
  basename(sub("^--file=", "", file_args[1])) == "run_morphometrics.R"
}

if (is_executed_as_script()) {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  run_analysis(opts)
}


