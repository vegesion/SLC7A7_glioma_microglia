# SLC7A7 in glioma-associated microglia

Analysis code for:

> **Tissue-level and cell-intrinsic inflammatory signals diverge with SLC7A7 expression in glioma-associated microglia**

This repository contains every script used to produce the figures and the
additional files of the manuscript. It contains **no data**. All primary data
are public; `docs/DATA_ACCESS.md` gives the accessions, versions and access
dates, and `docs/PIPELINE.md` gives the order in which the scripts run.

---

## What the analysis does

Four layers are analysed, and the unit of analysis is the donor or the tumor in
every one of them:

| Layer | Source | Unit |
|---|---|---|
| Bulk survival | TCGA-LGG, TCGA-GBM, CGGA-693, CGGA-325 (1,906 patients) | patient |
| Single cell | GBmap IDH-wildtype atlas (338,564 cells, 110 donors, 16 datasets) | donor |
| Spatial | GBM-space, 97 Visium sections from 12 tumors | tumor |
| Cell experiment | BV2 murine microglia, four-arm LPS design | biological replicate |

Three design choices run through the whole pipeline and explain most of what
the code does:

1. **Donor is the unit of analysis.** Cell-level tests are not used for
   inference anywhere.
2. **Expression is treated as continuous.** Cells are never dichotomised into
   positive and negative, because detection is confounded with sequencing
   depth.
3. **Sequencing depth is a competing hypothesis, not a nuisance.** Every
   donor-wise pathway analysis is repeated with depth as the predictor, and a
   result is carried forward only where the depth-only analysis fails to
   reproduce its direction.

---

## Running it

```r
# 1. restore the recorded environment (R 4.4.0, Bioconductor 3.20)
install.packages("renv")
renv::restore()                      # see docs/RENV_SETUP.md

# 2. point the project at your local copy of the primary data
file.copy(".Renviron.example", ".Renviron")   # then edit it

# 3. run the pipeline in order
source("analysis/01a_tcga_download_prepare.R")
# ... see docs/PIPELINE.md
```

Paths are resolved from the project root with `here::here()` and are set in
`config/config.R`. Large primary files (the GBmap `.h5ad`, the Visium objects,
the normalised count matrices) are **not** copied into the repository; the
location of your copy is read from `.Renviron`.

A single random seed, 42, is used throughout and is set in `config/config.R`.

---

## Layout

```
R/           function library, sourced by the analysis scripts
analysis/    the pipeline, numbered in execution order
figures/     scripts that only draw; they read results, they do not compute
config/      paths, parameters and thresholds -- the only place to edit them
docs/        pipeline, data access, environment, and the output map
data/        empty; primary data are not distributed here
```

`docs/OUTPUT_MAP.md` maps every additional file and every main figure to the
script that produces it. That is the fastest way to check a specific number in
the manuscript against the code.

---

## Analyses that were run and are not reported

Three analyses were run during the work and are deliberately absent from the
manuscript. They are named here so that their absence is not mistaken for an
oversight:

- **Mendelian randomisation.** Instrument strength was adequate (F > 89) but
  every estimate was null; the analysis cannot support a causal claim and is
  not reported.
- **Cell-cell communication (CellChat).** Its permutation test is at the level
  of single cells, which conflicts with design choice 1 above.
- **Cell-level differential expression on a positive/negative split.** The
  result was dominated by mitochondrial and ribosomal genes, the signature of a
  depth artefact, which is what motivated design choice 2.

A trajectory screen (`09_trajectory_exploratory.R`) is reported in the
manuscript as a descriptive figure only. No statistical test is applied to it
and no claim rests on it.

---

## Licence

Code is released under the MIT licence (`LICENSE`). The primary data remain
under the terms of their original repositories.
