Final Project – Automated Geometric Morphometric Analysis of Fish Vertebrae
John — MARB 6360


 Project Overview

This project automates a complete geometric morphometric workflow for quantifying vertebra shape differences between two fish populations (Matnog and Sacol). The analysis imports landmark coordinate files, performs Procrustes alignment, calculates PCA, compares mean shapes, evaluates allometry, and statistically tests for population-level morphological divergence.

The goal is to create a reproducible, fully automated pipeline for morphometric analysis while producing interpretable biological conclusions.









 Problem Statement
 
The objective is to determine whether two geographically distinct populations of fishes display:

Significant differences in vertebra shape

Differences in overall vertebra size

Shape differences that remain after removing size effects (allometry)

Because manually processing many vertebra images is inefficient and prone to error, an automated solution is needed.






Strategy Used

The workflow follows standard procedures in geometric morphometrics:

1. Load and assemble TPS landmark data

All *Back.TPS files from two population folders are discovered and read automatically.

Data are assembled into a unified Procrustes-friendly array.

2. Perform Generalized Procrustes Analysis (GPA)

GPA removes:

translation

rotation

scale

This isolates pure shape variation.

3. PCA visualization

Principal Component Analysis on Procrustes-aligned coordinates allows:

visualization of major axes of shape variation

detection of clustering by population

confidence ellipses

labeling specimens by file ID

4. Compare raw mean shapes

Group mean landmark configurations are computed and plotted together.

5. Quantify within-group variation

Procrustes distances to group means are used to compare morphological spread in each population.

6. Correct for allometry

A regression model (shape ~ centroid size) is used to compute size-corrected residual shape coordinates.

7. Re-align size-corrected means

To prevent rotation/mirroring artifacts, mean residual shapes are re-aligned using GPA.

8. Statistical testing

Three permutational Procrustes ANOVAs evaluate:

Test	Biological Question	Result
shape ~ group	Do populations differ in shape?	Yes (p ≈ 0.022)
size ~ group	Do populations differ in centroid size?	Yes (p ≈ 0.009)
shape ~ group + size	Does shape differ after size correction?	Yes (group p ≈ 0.011)
📊 Workflow Diagram
 ┌───────────────────────┐
 │  Detect TPS files      │
 │  in two populations    │
 └─────────┬─────────────┘
           │
           ▼
 ┌───────────────────────┐
 │ Read TPS landmarks     │
 │ Build 3D data array    │
 └─────────┬─────────────┘
           │
           ▼
 ┌────────────────────────┐
 │ GPA                    │
 │ Normalize shape data    │
 └─────────┬──────────────┘
           │
           ▼
 ┌────────────────────────┐
 │ PCA visualization       │
 │ + ellipses + labels     │
 └─────────┬──────────────┘
           │
           ▼
 ┌────────────────────────┐
 │ Compute raw mean shapes│
 │ Overlay Matnog/Sacol   │
 └─────────┬──────────────┘
           │
           ▼
 ┌────────────────────────────┐
 │ Procrustes distance metrics│
 │ Within-group variation      │
 └─────────┬──────────────────┘
           │
           ▼
 ┌────────────────────────┐
 │ Allometric regression   │
 │ shape ~ size            │
 └─────────┬──────────────┘
           │
           ▼
 ┌─────────────────────────────┐
 │ Compute size-corrected means │
 │ Re-align via GPA             │
 └─────────┬────────────────────┘
           │
           ▼
 ┌─────────────────────────────┐
 │ Statistical summary          │
 │ (3 Procrustes ANOVAs)        │
 └─────────────────────────────┘


📈 Summary of Key Results
1️⃣ Shape Differences Between Populations

p ≈ 0.022
✔ Significant shape divergence
✔ Population explains ~31% of total shape variation (R² ≈ 0.314)

Matnog and Sacol differ meaningfully in vertebra morphology.

2️⃣ Size Differences Between Populations

p ≈ 0.009
✔ Highly significant size difference
✔ Population explains ~84% of centroid size variation (R² ≈ 0.84)

Fish from the two populations differ substantially in vertebra size.

3️⃣ Shape Differences After Accounting for Size

group effect: p ≈ 0.011 → still significant

size effect: p ≈ 0.070 → marginal, not significant at α = 0.05

Shape differences persist even after removing size effects.
This demonstrates true morphological divergence rather than simple scaling differences.







Overall Biological Interpretation

The two populations exhibit significant vertebra shape differences.

They also differ strongly in size, but size alone does not explain the shape divergence.

After correcting for allometry, vertebra shape remains significantly different.

Conclusion:
The Matnog and Sacol fish populations show real morphological differentiation, consistent with localized adaptation, population structure, or ecological divergence.