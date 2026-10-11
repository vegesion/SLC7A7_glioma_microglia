# data/

This directory is intentionally empty. No primary data are distributed with
this repository.

```
data/
  raw/        primary downloads (TCGA, CGGA, GBmap h5ad, Visium objects)
  processed/  intermediate matrices and tables written by the pipeline
```

Both subdirectories are git-ignored. `docs/DATA_ACCESS.md` lists every source
with its accession, version and access date, and the scripts that populate
these directories are `01a`, `01b` and `06`.

If your copy of the primary data lives outside the project, set the search root
in `.Renviron` rather than editing paths inside the scripts; `config/config.R`
reads it from there.
