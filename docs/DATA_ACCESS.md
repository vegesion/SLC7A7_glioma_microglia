# Data access

No primary data are distributed with this repository. Every source below is
public. Accession, version and access date are given so that the exact
snapshot can be recovered.

## Bulk transcriptomes and clinical tables

| Cohort | Source | Version / access |
|---|---|---|
| TCGA-LGG | Genomic Data Commons, via `TCGAbiolinks` 2.34.1 | downloaded 25 July 2026 |
| TCGA-GBM | Genomic Data Commons, via `TCGAbiolinks` 2.34.1 | downloaded 25 July 2026 |
| CGGA-693 | Chinese Glioma Genome Atlas, `mRNAseq_693` (batch 1) | downloaded 28 July 2026 |
| CGGA-325 | Chinese Glioma Genome Atlas, `mRNAseq_325` (batch 2) | downloaded 28 July 2026 |

`01a_tcga_download_prepare.R` performs the GDC download; see
`docs/GDC_DOWNLOAD.md` for the client setup. The CGGA files are downloaded by
hand from the CGGA portal and placed in `data/raw/`.

## Single-cell atlas

GBmap, the IDH-wildtype core atlas.

- CELLxGENE collection `999f2a15-3d7e-440b-96ae-2c806799c08c`
- dataset `c888b684-6c51-431f-972a-6c963044cef0`
- version `861acfd8-25f0-418b-a445-aa96da232827`
- downloaded 25 July 2026

The `.h5ad` is read directly with `hdf5r`, without `reticulate`, by
`06_sc_preprocess_gbmap.R`. The expression matrix distributed under the
CELLxGENE schema is **not** batch-corrected; the integration applies to the
embedding only. Every expression-based analysis therefore carries a dataset
random effect.

## Spatial transcriptomes

GBM-space. Processed Visium objects downloaded from www.gbmspace.org on
26 July 2026. The underlying series is deposited in ArrayExpress under
accession **E-MTAB-17164**. Ninety-seven sections from twelve tumors.

## Where to put the files

Set the search root in `.Renviron` (copy `.Renviron.example`). `config/config.R`
reads it and resolves every path from there, so no script contains an absolute
path. If the files are inside the project, `data/raw/` and `data/processed/`
are searched first and nothing needs to be set.
