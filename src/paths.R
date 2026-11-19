# paths.R : where files live, on Lifebit.
# The whole folder sits at /home/jovyan/session_data/workflow-007 and everything is relative to it.
# Public-data notebooks load it with:  source("/home/jovyan/session_data/workflow-007/src/paths.R")

ROOT <- "/home/jovyan/session_data/workflow-007"   # the repo folder itself
DATA <- file.path(ROOT, "data")                  # downloaded public files:  data/<study>/<file>
NCBI <- file.path(DATA, "ncbi")                  # NCBI gene table shared by the studies
GENES <- file.path(ROOT, "genes")                # gene-list files (one symbol per line)
RESULTS <- file.path(ROOT, "results")            # everything the NIH notebooks write (results/rnaseq_out_v2/...)
SEED <- 20260928                                 # random seed for any method that draws random numbers

.in <- function(base, study, file) {             # make the folder if needed, return folder or file path
  d <- file.path(base, study)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  if (is.null(file)) d else file.path(d, file)
}
raw <- function(gse, file = NULL) .in(DATA, gse, file)                              # raw("GSE65391", "x.gz") -> .../data/GSE65391/x.gz
art <- function(gse, file = NULL) .in(file.path(DATA, "run_artifacts"), gse, file)  # results:  .../data/run_artifacts/GSE65391/x.rds
fig <- function(gse, file = NULL) .in(file.path(ROOT, "figures"), gse, file)        # figures:  .../figures/GSE65391/x.png

read_gene_set <- function(name) {                # lines starting with # are notes
  x <- trimws(readLines(file.path(GENES, name)))
  x[x != "" & !startsWith(x, "#")]
}

dir.create(NCBI, recursive = TRUE, showWarnings = FALSE)
