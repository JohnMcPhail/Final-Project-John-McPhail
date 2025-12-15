# Automated Geometric Morphometric Analysis of Fish Vertebrae

This repository contains an automated geometric morphometrics workflow to
compare vertebra shape between two fish populations (**Matnog** vs **Sacol**)
using landmark coordinates stored in TPS files.

## Repo Organization

- `data/raw/` — raw input data (TPS landmark files and associated X-ray images)
- `data/processed/` — derived/processed outputs written by the pipeline
- `scripts/` — runnable analysis scripts
- `results/` — analysis outputs (summary + figures)
- `reports/` — optional R Markdown report

## Data Notes

The analysis reads `*_Back.TPS` landmark files from:

- `data/raw/johns_landmarked_photos/2024-12_nmnh_c-matnog_a-matnog-sorsogon_xrays/xrays_cropped-for-geomorph/`
- `data/raw/johns_landmarked_photos/2024-12_nmnh_c-taluksangay_a-sacol-island_xrays/xrays_cropped-for-geomorph/`

The TPS files are the **raw analysis inputs** for this pipeline (landmarks were
created upstream in a separate landmarking step).

## How To Run

From the repository root:

```bash
Rscript scripts/run_morphometrics.R
```

To see options:

```bash
Rscript scripts/run_morphometrics.R --help
```

### Requirements

- R (tested with R 4.5.x)
- R package: `geomorph` (plus its dependencies)

## Outputs

After running `scripts/run_morphometrics.R`, the pipeline writes:

- `results/summary.txt` — model outputs (shape ~ group, size ~ group, shape ~ group + size)
- `results/figures/` — automated visualizations (PCA, mean shapes, variation boxplot)
- `data/processed/specimens.csv` — specimen IDs, groups, TPS paths
- `data/processed/geomorph_data.rds` — saved R objects (coords, centroid sizes, PCA scores)

## Report (Optional)

`reports/morphometrics_report.Rmd` can be knitted after running the pipeline to
embed the generated figures and summary into an HTML report.
