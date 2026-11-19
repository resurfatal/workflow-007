# ===== CELL 1 =====
library(readxl)

# ---- settings (change here, nowhere else) -----------------------------------
OUT        <- "/mnt/filesystems/nanostring_out"            # where notebook 01 wrote
MASTER_CSV <- file.path(OUT, "nanostring_master.csv")
JDM_REGEX  <- "^JDM IFN scores OADR.*\\.xlsx$"             # searched for below
JDM_SEARCH <- c("/mnt/filesystems", "/mnt/session_data", "/home/jovyan")
JDM_XLSX   <- NULL        # or a full path, to skip the search

# Which sheet is which version. Leave empty and the notebook guesses from the sheet
# name (raw / STD / "no norm"). If it guesses wrong, fill in, e.g.
#   SHEET_VERSION <- c("Sheet1" = "raw", "Sheet2" = "STD", "Sheet3" = "raw_nonorm")
SHEET_VERSION   <- c()

ID_COL_OVERRIDE <- "study sample"   # JDM codes like 015-01-01 = patient 015, visit 01, sample 01 (set NULL to auto-guess)
ADD_OWN_HK_NORM <- TRUE   # Part 6

PANEL <- c("CXCL10","DDX60","EPSTI1","GBP1","HERC5","HERC6","IFI27","IFI44","IFI44L",
           "IFI6","IFIT1","IFIT2","IFIT3","IFIT5","IFNA2","IFNB1","IRF7","ISG15",
           "LAMP3","LY6E","MX1","OAS1","OAS2","OAS3","OASL","RSAD2","RTP4","SAMD9",
           "SIGLEC1","SOCS1","SPATS2L","USP18")
HOUSEKEEPING <- c("ALAS1","HPRT1","TBP","TUBB")

normalize_id <- function(x) {                      # same rule as notebook 01: O-6033 -> O6033, 1013.0 -> 1013
  s <- trimws(as.character(x))
  s[is.na(x)] <- ""
  s[tolower(s) %in% c("nan", "none", "nat", "na")] <- ""
  s <- sub("\\.0$", "", s)
  toupper(gsub("[-_ ]", "", s))
}
nc <- function(x) gsub("[^a-z0-9]", "", tolower(x))
GENE_BY_NC <- setNames(c(PANEL, HOUSEKEEPING), nc(c(PANEL, HOUSEKEEPING)))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# rbind that tolerates different columns (missing ones become NA)
bind_fill <- function(dfs) {
  cols <- unique(unlist(lapply(dfs, names)))
  do.call(rbind, lapply(dfs, function(d) {
    for (m in setdiff(cols, names(d))) d[[m]] <- NA
    d[, cols, drop = FALSE]
  }))
}

# ===== CELL 2 =====
master <- read.csv(MASTER_CSV, check.names = FALSE, stringsAsFactors = FALSE)
master$dataset      <- "RGM"
master$version      <- "raw_counts"
master$sample_label <- if ("raw_filename" %in% names(master)) master$raw_filename else ""
master$patient_key  <- if ("patient_id" %in% names(master)) normalize_id(master$patient_id) else ""

cat("rows x columns:", nrow(master), "x", ncol(master), "\n")
if ("record_type" %in% names(master)) print(table(master$record_type))
cat("panel + housekeeping genes present:", sum(c(PANEL, HOUSEKEEPING) %in% names(master)),
    "of", length(c(PANEL, HOUSEKEEPING)), "\n")
head(master[, intersect(c("patient_key","diagnosis1","specimen_type","record_type","qc_pass","run_date"), names(master))], 5)

# ===== CELL 3 =====
if (is.null(JDM_XLSX)) {
  hits <- unlist(lapply(JDM_SEARCH[dir.exists(JDM_SEARCH)], function(d)
    list.files(d, pattern = JDM_REGEX, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)))
  hits <- hits[!grepl("~\\$", hits)]
  if (length(hits) == 0) stop("No file like 'JDM IFN scores OADR.xlsx' under ", paste(JDM_SEARCH, collapse = ", "),
                              ". Set JDM_XLSX in Part 1.")
  JDM_XLSX <- hits[1]
}
cat("using:", JDM_XLSX, "\n")
sheet_names <- excel_sheets(JDM_XLSX)
cat("sheets:", paste(sheet_names, collapse = " | "), "\n")

ID_WORDS <- c("ptid","patientid","enrollid","nihenrollid","nihenrollmentid","enrollmentid","id","sampleid",
              "subjectid","subject","patient","ptidnihenrollid","studyid","mrn")

# 'IFI27 STD' -> list(gene="IFI27", suffix="STD");  'IFI27' -> suffix "";  anything else -> NULL
# Only a known version word may follow the gene name. So 'CXCL10/IP-10' (some other
# measurement) is NOT treated as the CXCL10 gene column.
SUFFIX_OK <- "^(raw|std|standard|norm|normalized|normalised|nonorm|nonnorm|unnorm|log2|log)*$"
split_gene_col <- function(x) {
  t <- trimws(as.character(x))
  if (is.na(t) || t == "") return(NULL)
  m <- regmatches(t, regexec("^([A-Za-z0-9]+)[ _.-]*(.*)$", t))[[1]]
  if (length(m) == 0) return(NULL)
  k <- nc(m[2])
  if (!(k %in% names(GENE_BY_NC))) return(NULL)
  if (!grepl(SUFFIX_OK, nc(m[3]))) return(NULL)
  list(gene = unname(GENE_BY_NC[k]), suffix = trimws(m[3]))
}

guess_header_row <- function(raw) {
  best <- 1; best_n <- -1
  for (i in seq_len(min(15, nrow(raw)))) {
    vals <- unlist(raw[i, ], use.names = FALSE)
    n <- sum(vapply(vals, function(v) !is.null(split_gene_col(v)), logical(1)))
    if (n == 0) n <- 0.01 * sum(!is.na(vals))
    if (n > best_n) { best <- i; best_n <- n }
  }
  best
}

guess_version <- function(sheet) {
  s <- tolower(sheet); tags <- character(0)
  if (grepl("raw", s)) tags <- c(tags, "raw")
  if (grepl("std|standard", s)) tags <- c(tags, "STD")
  if (grepl("no[ _-]*norm|non[ _-]*norm|un[ _-]*norm|no[ _-]*form|not[ _-]*norm", s)) tags <- c(tags, "nonorm")
  if (length(tags)) paste(tags, collapse = "_") else gsub("^_|_$", "", gsub("[^A-Za-z0-9]+", "_", sheet))
}

jdm_sheets <- list()
for (sh in sheet_names) {
  raw <- as.data.frame(suppressMessages(read_excel(JDM_XLSX, sheet = sh, col_names = FALSE, col_types = "text")),
                       stringsAsFactors = FALSE)
  raw <- raw[rowSums(!is.na(raw)) > 0, , drop = FALSE]
  if (nrow(raw) < 2) { cat("\n[", sh, "] empty, skipped\n", sep = ""); next }
  h <- guess_header_row(raw)
  hdr <- unlist(raw[h, ], use.names = FALSE)
  hdr <- ifelse(is.na(hdr), paste0("unnamed_", seq_along(hdr)), trimws(hdr))
  df <- raw[-seq_len(h), , drop = FALSE]
  names(df) <- make.unique(hdr)
  rownames(df) <- NULL

  parsed    <- lapply(names(df), split_gene_col); names(parsed) <- names(df)
  gene_cols <- parsed[!vapply(parsed, is.null, logical(1))]
  # Which column holds the patient ID? Test every non-gene column against the patient IDs
  # of your RGM set (after the same cleaning) and take the one that matches most.
  rcc_keys <- unique(master$patient_key[master$patient_key != ""])
  cands <- setdiff(names(df), names(gene_cols))
  cands <- cands[!grepl("^unnamed_", cands)]
  cands <- cands[!grepl("ifn|score|geomean|index|std|norm|gene$", cands, ignore.case = TRUE) | grepl("study id|nano id", cands, ignore.case = TRUE)]  # score columns are never IDs
  ov <- vapply(cands, function(cn) {
    k <- unique(normalize_id(df[[cn]])); sum(k[k != ""] %in% rcc_keys) }, numeric(1))
  idc <- if (!is.null(ID_COL_OVERRIDE) && ID_COL_OVERRIDE %in% names(df)) ID_COL_OVERRIDE else if (length(ov) && max(ov) > 0) {
    names(ov)[which.max(ov)] } else {
    w <- cands[nc(cands) %in% c("studyid","nanoid","ptid","patientid","enrollid","nihenrollid","enrollmentid","sampleid")]
    if (length(w)) w[1] else cands[1] }
  top <- sort(ov, decreasing = TRUE); top <- top[top > 0]
  cat(sprintf("\n[%s] patient-ID candidates (columns whose values match your RGM patient IDs):\n", sh))
  if (length(top)) print(head(top, 5)) else cat("   none match - falling back to a column named like an ID; set ID_COL_OVERRIDE if wrong\n")
  other     <- setdiff(names(df), c(names(gene_cols), idc))
  other     <- other[!grepl("^unnamed_", other)]
  score_cols <- other[grepl("ifn|score|geomean|index|std|norm", other, ignore.case = TRUE) & !grepl("id$", other, ignore.case = TRUE)]
  ver <- if (sh %in% names(SHEET_VERSION)) unname(SHEET_VERSION[sh]) else guess_version(sh)

  jdm_sheets[[sh]] <- list(df = df, id = idc, genes = gene_cols, scores = score_cols, version = ver)
  cat(sprintf("\n[%s] %d rows x %d cols | header row %d | version = '%s'\n", sh, nrow(df), ncol(df), h, ver))
  cat("   ID column chosen:", idc, " e.g.", paste(head(df[[idc]], 3), collapse = ", "),
      sprintf(" (%d of %d values match RGM patients)\n", sum(normalize_id(df[[idc]]) %in% rcc_keys), nrow(df)))
  cat("   gene columns:", length(gene_cols), " e.g.", paste(head(names(gene_cols), 4), collapse = ", "), "\n")
  cat("   score columns carried along:", if (length(score_cols)) paste(score_cols, collapse = ", ") else "none", "\n")
  cat("   other columns left behind:", paste(head(setdiff(other, score_cols), 12), collapse = ", "), "\n")
}

# ===== CELL 4 =====
# JDM codes look like "015-01-01" (patient-visit-sample) or just "103". The key is "JDM" + the patient
# part ("JDM015"). The prefix is deliberate: JDM study codes are a different numbering from NIH
# enrollment IDs, so they must never match an RGM ID by accident (e.g. JDM "103" vs RGM "103").
jdm_code  <- function(x) { x <- trimws(as.character(x)); x[is.na(x)] <- ""; sub("-.*$", "", x) }
jdm_key   <- function(x) { c <- jdm_code(x); ifelse(c == "", "", paste0("JDM", c)) }
jdm_visit <- function(x) { x <- trimws(as.character(x)); ifelse(grepl("-", x), sub("^[^-]*-([^-]*).*$", "\\1", x), NA) }

jdm_blocks <- list()
for (sh in names(jdm_sheets)) {
  d <- jdm_sheets[[sh]]; df <- d$df
  # one block per distinct suffix ('IFI27 raw' and 'IFI27 STD' in one sheet become two versions)
  suffixes <- unique(vapply(d$genes, function(g) g$suffix, character(1)))
  for (suf in suffixes) {
    cols <- names(d$genes)[vapply(d$genes, function(g) g$suffix == suf, logical(1))]
    blk <- data.frame(dataset = "JDM",
                      version = if (suf == "") d$version else paste(d$version, suf, sep = "_"),
                      patient_key  = jdm_key(df[[d$id]]),
                      sample_label = as.character(df[[d$id]]),
                      stringsAsFactors = FALSE)
    for (cc in cols) blk[[d$genes[[cc]]$gene]] <- suppressWarnings(as.numeric(df[[cc]]))
    blk$jdm_sheet <- sh
    blk$jdm_visit <- jdm_visit(df[[d$id]])
    # the workbook stores dates as Excel serial numbers (37839 = 2003-08-06)
    if ("sample date" %in% names(df))
      blk$sample_date <- format(as.Date(suppressWarnings(as.numeric(df[["sample date"]])), origin = "1899-12-30"))
    for (sc in d$scores) blk[[paste0("jdm_", gsub("^_|_$", "", gsub("[^A-Za-z0-9]+", "_", sc)))]] <- df[[sc]]
    jdm_blocks[[length(jdm_blocks) + 1]] <- blk
  }
}
jdm_long <- bind_fill(jdm_blocks)
cat("JDM rows:", nrow(jdm_long), "\n"); print(table(jdm_long$version))

combined <- bind_fill(list(master, jdm_long))
front <- c("dataset", "version", "patient_key", "sample_label")
combined <- combined[, c(front, setdiff(names(combined), front))]

stopifnot("a birth-date column reached the output" = !any(grepl("birth|\\bdob\\b", names(combined), ignore.case = TRUE, perl = TRUE)))
out_path <- file.path(OUT, "nanostring_all_datasets_long.csv")
write.csv(combined, out_path, row.names = FALSE, na = "")
cat(sprintf("saved %s: %d rows x %d columns\n", out_path, nrow(combined), ncol(combined)))
as.data.frame(table(dataset = combined$dataset, version = combined$version))

# ===== CELL 5 =====
shared_genes <- intersect(c(PANEL, HOUSEKEEPING), intersect(names(master), names(jdm_long)))
miss_in_jdm  <- setdiff(intersect(c(PANEL, HOUSEKEEPING), names(master)), names(jdm_long))
cat("5a. genes in both:", length(shared_genes), "| in RGM set but not in JDM workbook:",
    if (length(miss_in_jdm)) paste(miss_in_jdm, collapse = ", ") else "none", "\n")

med_log2 <- function(v) median(log2(pmax(as.numeric(v), 0) + 1), na.rm = TRUE)
# like-for-like comparison: whole-blood patient and healthy-control samples only
# (the 7,690 files also hold standards, in-vitro experiments, BAL/PBMC, covid, test runs...)
rcc <- master
if ("record_type" %in% names(master))  rcc <- rcc[rcc$record_type %in% c("patient_sample", "healthy_control"), ]
if ("specimen_type" %in% names(master)) rcc <- rcc[rcc$specimen_type %in% "whole_blood", ]
cat("RGM rows used for 5b:", nrow(rcc), "of", nrow(master), "\n")
versions <- sort(unique(jdm_long$version))
scale_tab <- data.frame(gene = shared_genes, RGM_raw = sapply(shared_genes, function(g) med_log2(rcc[[g]])))
for (v in versions) {
  sub <- jdm_long[jdm_long$version == v, ]
  scale_tab[[paste0("JDM_", v)]] <- sapply(shared_genes, function(g) med_log2(sub[[g]]))
}
rownames(scale_tab) <- NULL
cat("\n5b. median log2(count+1) per gene\n")
print(head(format(scale_tab, digits = 3), 40))

jcols <- grep("^JDM_", names(scale_tab), value = TRUE)
if (length(jcols)) {
  par(mfrow = c(1, length(jcols)), mar = c(4, 4, 2, 1))
  lim <- c(0, max(unlist(scale_tab[, -1]), na.rm = TRUE) + 1)
  for (v in jcols) {
    plot(scale_tab$RGM_raw, scale_tab[[v]], xlim = lim, ylim = lim, pch = 19, cex = .8,
         xlab = "RGM set, raw (log2)", ylab = paste(v, "(log2)"), main = v)
    abline(0, 1, lty = 2)
  }
  par(mfrow = c(1, 1))
}

mk   <- unique(master$patient_key[master$patient_key != ""])
jk   <- unique(jdm_long$patient_key[jdm_long$patient_key != ""])
cat(sprintf("\n5c. patients: RGM set %d | JDM %d | ID overlap %d (JDM codes are a separate numbering, so ~0 is expected)\n",
            length(mk), length(jk), length(intersect(mk, jk))))

# The real test of "same samples": do the collection dates even overlap?
rng <- function(x) { d <- as.Date(x[!is.na(x) & x != ""]); if (length(d)) range(d) else as.Date(c(NA, NA)) }
r_rcc <- rng(master$sample_date); r_jdm <- rng(jdm_long$sample_date)
cat("   RGM set sample dates:", format(r_rcc[1]), "to", format(r_rcc[2]), "\n")
cat("   JDM sample dates:    ", format(r_jdm[1]), "to", format(r_jdm[2]), "\n")
if ("run_date" %in% names(master)) { rr <- rng(master$run_date); cat("   RGM set RUN dates:   ", format(rr[1]), "to", format(rr[2]), "\n") }
if (!is.na(r_rcc[1]) && !is.na(r_jdm[2]) && r_jdm[2] < r_rcc[1]) {
  cat("   JDM samples were all collected BEFORE any RGM-set sample: no sample can be in both sets,\n",
      "  so there are no technical replicates and RUV-III has nothing to anchor on.\n")
}

# ===== CELL 6 =====
if (ADD_OWN_HK_NORM) {
  if (all(HOUSEKEEPING %in% names(combined))) {
    raw_like <- !grepl("STD", combined$version, ignore.case = TRUE)
    hkm <- sapply(HOUSEKEEPING, function(g) suppressWarnings(as.numeric(combined[[g]])))
    ok  <- raw_like & complete.cases(hkm) & apply(hkm, 1, function(r) all(r > 0))
    geo <- exp(rowMeans(log(hkm)))
    scale_f <- mean(geo[ok])
    for (g in setdiff(shared_genes, HOUSEKEEPING)) {
      v <- suppressWarnings(as.numeric(combined[[g]]))
      col <- rep(NA_real_, nrow(combined))
      col[ok] <- log2(pmax(v[ok] / geo[ok] * scale_f, 0) + 1)
      combined[[paste0(g, "_hk")]] <- col
    }
    cat(sprintf("added _hk columns for %d rows; skipped %d raw-like rows with a zero/missing housekeeping gene and %d STD rows\n",
                sum(ok), sum(raw_like & !ok), sum(!raw_like)))
    write.csv(combined, out_path, row.names = FALSE, na = "")
    cat("re-saved", out_path, "\n")
  } else {
    cat("housekeeping genes missing from the table:", paste(setdiff(HOUSEKEEPING, names(combined)), collapse = ", "), "\n")
  }
}
