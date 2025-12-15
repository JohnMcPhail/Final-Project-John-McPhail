# MODIFICATIONS (by Codex / GPT-5.2)

I am **Codex using GPT-5.2**. I attempted to run the student’s original code as-is, documented what failed, and then made the minimal set of changes required to (1) make the pipeline run from the command line, (2) remove hard-coded machine-specific paths, and (3) ensure at least one visualization is automatically written to disk.

## What Failed Initially (Reproduction)

Running the original script failed outside of RStudio:

```bash
Rscript final_coding_project_computational_biology.R
```

Error observed:

```text
Error: RStudio not running
Execution halted
```

Root cause: the script depended on `rstudioapi::getActiveDocumentContext()` to set the working directory.

## Repo Organization Changes

### 1) Created a standard project layout

Created directories to match common “raw → processed → results” organization:

- `data/raw/`
- `data/processed/`
- `scripts/`
- `results/figures/`
- `reports/`

Added placeholder files so empty directories can exist in git:

- `data/processed/.gitkeep`
- `results/.gitkeep`
- `results/figures/.gitkeep`

Why: makes the project easier to understand and run; matches typical computational biology repo structure.

### 2) Moved/renamed the raw data folder

Moved:

- `Johns landmarked photos/` → `data/raw/johns_landmarked_photos/`

Why: keeps raw inputs under `data/raw/` and removes a space-heavy path from the analysis defaults.

## Script and Report Changes

### 3) Renamed and relocated the main script

Moved/renamed:

- `final_coding_project_computational_biology.R` → `scripts/run_morphometrics.R`

Why: puts runnable code in `scripts/` and gives it an “entrypoint” name that matches what it does.

### 4) Removed RStudio-only dependency and hard-coded working-directory logic

Removed:

```r
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
```

Why: `rstudioapi` fails in non-RStudio contexts, breaking command-line reproducibility.

### 5) Removed unused heavy dependency

Removed:

```r
library(tidyverse)
```

Why: the pipeline does not use tidyverse functions; removing it reduces install burden and avoids namespace conflicts.

### 6) Parameterized inputs/outputs and added `--help`

Added a simple CLI argument parser so the pipeline can be run from repo root and optionally pointed at other folders:

```r
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
  ...
}
```

Why: avoids editing code to run the analysis and supports reproducible “one command” execution.

### 7) Fixed `procD.lm` evaluation inside a function

Inside `run_analysis()`, `geomorph::procD.lm()` failed with:

```text
Error: Cannot find data in global environment.
```

Fix: provide an explicit `data=` list to `procD.lm` so the model can find `coords`, `CS`, and `group` when run inside a function:

```r
model_data <- list(coords = coords, CS = CS, group = group)
shape_group_fit <- geomorph::procD.lm(coords ~ group, data = model_data, iter = opts$iter, print.progress = FALSE)
```

Why: prevents scoping issues and makes the script reliable when wrapped in functions (and when run via `Rscript`).

### 8) Replaced hard-coded p-values/results with computed results

The original script printed fixed “p ≈ ...” values via `cat(...)`.

Fix: compute the stats directly from model output and write them to `results/summary.txt`:

```r
shape_group_anova <- anova(shape_group_fit)
sprintf("   - p = %s", fmt_p(get_anova_value(shape_group_anova, "group", "Pr(>F)")))
```

Why: prevents reports from drifting out of sync with the actual data/models being run.

### 9) Automated figure creation (writes files)

Added a helper and used it to save plots as PNGs:

```r
save_png(file.path(fig_dir, "pca_pc1_pc2.png"), plot_fun = function() { ... })
```

Why: fulfills the “automate at least one visualization” requirement in a way that works non-interactively.

### 10) Wrote processed outputs for reproducibility

Added:

- `data/processed/specimens.csv` (metadata table)
- `data/processed/geomorph_data.rds` (saved R objects)

Why: preserves the “wrangled” dataset that downstream steps can reuse without re-reading all TPS files.

### 11) Replaced the Rmd with a minimal, portable report

Moved/renamed:

- `finalproject.Rmd` → `reports/morphometrics_report.Rmd`

Replaced machine-specific Windows paths with instructions and optional embedding of generated outputs.

Why: the original Rmd hard-coded absolute Windows paths and duplicated the whole script; the updated report points to the script and reads outputs from `results/`.

## Documentation and Git Hygiene

### 12) Updated `.gitignore`

Updated `.gitignore` to:

- ignore Windows `*:Zone.Identifier`
- ignore R session artifacts (`.RData`, `.Rhistory`)
- ignore generated outputs under `results/` and `data/processed/` while keeping `.gitkeep`
- support both old and new raw-data folder locations

Why: keeps the repo clean and prevents accidentally committing generated artifacts or OS metadata.

### 13) Removed tracked OS/RStudio artifacts

Removed tracked files:

- `.RData`
- `.Rhistory`
- `README.md:Zone.Identifier`
- `final_coding_project_computational_biology.R:Zone.Identifier`
- `finalproject.Rmd:Zone.Identifier`

Why: these are not part of the analysis and commonly cause noise/merge conflicts.

### 14) Rewrote `README.md`

Replaced the original README with a concise, “how to run / where outputs go” README aligned to the new structure.

Why: graders and future users need quick, correct execution instructions and a clear directory map.

