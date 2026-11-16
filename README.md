# endotype-transcripts-007

Tanisha's copy of the method in Anne Deslattes Mays' repository (NIH-NLM/endotypes-transcriptomics), applied to the internal RNA-seq data (RGM: four sequencing facilities) and the JDM-lab data, run on Lifebit.

## Start here

## Setup: the terminal commands, in order

Run these in a Lifebit **terminal** (not a notebook cell). Steps 1 to 3 are done once per Lifebit session; the environment from step 2 stays in your home folder.

**1. Put the folder in place.** It must be at exactly `/home/jovyan/session_data/endotype-transcripts-007`. Every notebook reads and writes under that folder (`genes/`, `src/paths.R`, `data/` for public downloads, `results/` for outputs), so there is no copy step.

```bash
cd ~/session_data
ls endotype-transcripts-007            # should list: README.md OUTPUTS.md ipynb src genes env
```

If you already ran notebooks into the old `Attempt0` folder, move that work across once so nothing re-runs or is orphaned:

```bash
cd ~/session_data
[ -d Attempt0/endotype_results ] && mv Attempt0/endotype_results endotype-transcripts-007/results
[ -d Attempt0/data ] && mv Attempt0/data endotype-transcripts-007/data
```

**2. Build the environment** (R with limma, edgeR, sva, plus Python for notebooks 01 and 02). Takes several minutes; leave the terminal open until it prints "Done".

```bash
bash /home/jovyan/session_data/endotype-transcripts-007/env/setup_endotypes_env.sh
```

Every package line near the end should say `OK`. The environment is created at `/home/jovyan/envs/endotypes`. Do not use `conda install -p /home/jovyan/envs/wgcna ...`: that environment does not exist in the current Lifebit session.

**3. Check your input files are present.**

```bash
ls /mnt/filesystems | grep -i -E "Goldbach|HC_blood|MetadataPatients|4b_rpkm"
```

You need the Goldbach-Mansky workbook, `HC_blood_all_samples.xlsx`, `MetadataPatients_NoPHI.xlsx` and `4b_rpkm.csv`. (04 needs no JDM key sheet.)

**4. In the browser:** refresh the tab, then open each notebook and choose Kernel > Change Kernel:

| notebooks | kernel |
|---|---|
| `rnaseq/nih` 01, 02 | **Python (endotypes)** |
| `rnaseq/nih` 03 to 10 | **R (endotypes)** (the notebooks say "R (WGCNA)" in their text; ignore that) |
| `rnaseq/public` (all) | **R (endotypes)** |
| `nanostring` 01, 02 | **Python (endotypes)** |

Check R once, in a cell on the R kernel. All three must say `TRUE`:

```r
sapply(c("limma", "edgeR", "sva"), requireNamespace, quietly = TRUE)
```

**5. Run each notebook with Restart & Run All, in the order below.** Each reads the previous ones' files, so stop at the first error.

**Handy checks while you work** (terminal):

```bash
# what has been written so far
ls /home/jovyan/session_data/endotype-transcripts-007/results/rnaseq_out_v2/
# did 04 finish? (should list 9 files)
ls /home/jovyan/session_data/endotype-transcripts-007/results/rnaseq_out_v2/all_sources_correction/
# which kernels exist
jupyter kernelspec list
```

If something fails, copy the red error message and the cell it came from. For install problems, copy the last 30 lines of the terminal.

**Known harmless messages:** ComBat's "Found N batches / Standardizing Data" lines in 04; and in 02 the Checkpoint 10 histogram error is already fixed in this folder's copy.

Every file each notebook writes, and what it is for, is listed in **`OUTPUTS.md`**.

## Run order

| step | notebook (in `ipynb/rnaseq/nih/`) | language | what it does | writes (in `results/rnaseq_out_v2/`) |
|---|---|---|---|---|
| 01 | `01_facility_csvs` | Python | each facility tab to its own CSV, healthy controls added to the facility they were sequenced at | `facility_csv/` |
| 02 | `02_parse_metadata` | Python | whole blood, library size, patient metadata, drops post-treatment, one sample per patient | counts, sample sheets, `batch_by_diagnosis_eligible.csv` |
| 03 | `03_batch_combatseq` | R | the four facilities corrected with ComBat-seq on counts (RGM only) | `batch_correction/` |
| 04 | **`04_batch_all_sources`** | R | **RGM + JDM lab on one scale, with proof** | `all_sources_correction/` |
| 05 | `05_differential_expression` | R | each disease against healthy, limma-voom, raw counts, batch as a model term | `differential_expression/` |
| 06 | `06_clusters_heatmaps` | R | heatmaps and sample groups (RGM) | `clusters/` |
| 07 | `07_jdm_comparison` | R | JDM against each RGM disease, each against its own lab's healthy controls | `jdm_comparison/` |
| 08 | `08_sle_gse65391_vs_jdm` | R | SLE (public array GSE65391) against its own healthy children, then against every disease and JDM, with and without the interferon genes. Needs the three `public/GSE65391_*` notebooks first | `sle_comparison/` |
| 09 | `09_wgcna_modules` | R | gene modules (WGCNA), one network for RGM and one for the JDM lab, then module preservation both ways and module overlap | `wgcna/` |
| 10 | `10_sle_gse135779_vs_jdm` | R | the second SLE study (single-cell PBMC as pseudobulk) against its own healthy children, then against every disease, JDM and the first SLE study. Downloads 1.3 GB on the first run; needs the by-hand supplement | `sle2_comparison/` |

The numbers restart in each folder, so "notebook 05" in any text means `nih/05_differential_expression`. Notebooks in the other two folders are always named with their folder or study.

The two NanoString notebooks (`ipynb/nanostring/`: `01_rcc_parse_qc`, then `02_append_jdm`, which needs 01's `nanostring_master.csv`) and the notebooks in `ipynb/rnaseq/public/` are independent of this chain, except that 08 needs the three GSE65391 ones. Notebook 10 downloads and builds its own GSE135779 data, so the four `GSE135779_*` notebooks are Anne's separate version of the same steps and are not needed for it.

## How to read 04 in three questions

1. **Part 4: which sources can be corrected?** A source can be corrected when it holds healthy controls, or a disease, that another source also holds. That shared group is the ruler its shift is measured from. A source holding one diagnosis found nowhere else cannot be corrected by any method, because "this source" and "this disease" are the same thing. Such sources are left out of the combined table, listed with the reason, and never silently included.
2. **Part 6a: did the source pattern fall?** The share of the main patterns explained by *source* should drop from stage 1 to stage 3, while the share explained by *diagnosis* stays.
3. **Part 6b and 6c: is the biology still there?** The interferon genes must still be up in CANDLE, SAVI, AGS, UAID and JDM after correction (6b), and a disease found in two sources must give the same fingerprint in both (6c). The last cell prints a SUMMARY with PASS or CHECK for each.

## How this maps to Anne's repository

| here | Anne's notebook | note |
|---|---|---|
| 01, 02, 03 | 40 | reading the four facility files, ComBat-seq on counts |
| 05 | 41, 46 | limma-voom, patient as block, batch in the design |
| 06 | 42 | heatmaps of top genes, interferon genes, Ward trees |
| 07 | 46, 46b | each disease against its own healthy controls |
| **04** | **48, 49** | Anne tried to remove the lab difference with control genes (RUVg) and concluded that no setting works, because each lab held one kind of sample. 04 therefore anchors on the healthy controls that the sources really share, and reuses her **positive and negative controls** (28 interferon genes up, 11 NF-κB genes unchanged) as pass or fail checks. |
| 10 | 20 to 22, 48 | the second SLE study: download, metadata, pseudobulk, SLE vs healthy with batch in the model |
| 09 | 12 to 15, 25 to 26 | WGCNA per study and module preservation (the array modules are not projected yet) |
| 08 | 06c, 48 | SLE against its own healthy children, then compared with the other diseases |

## Results so far (what the data said)

These are readings of the printed output, from notebooks that ran on Lifebit. They are not final conclusions.

- **04, all-sources correction.** Six sources were corrected (JDM-lab, NCI-FNL, NIAID-NCI, NIAMS-Batch1, NIAMS-Batch2, NIAMS-unlabelled; 218 samples, 11,145 genes). **BCM-HGSC was left out**: it has only 8 samples (2 each of AGS, CANDLE, HC, SAVI), too few to measure its shift. Source explained 0.98 of PC1+PC2 before correction and 0.10 after. The interferon check passed for all 5 diseases (5 PASS, 0 FAIL) and the shared-disease check passed (5 PASS, 0 CHECK). The one CHECK was the crude 6a rule (diagnosis share fell 0.41 to 0.32 against a threshold of 0.33); source fell as intended and the other two checks passed, so it was treated as a false alarm.
- **07, JDM against each disease.** Rank correlation of fold changes: UAID 0.55, CANDLE 0.46, SAVI 0.35, AGS -0.32. Shared top genes all go in the same direction for CANDLE, SAVI and UAID. Interferon score (each against its own lab's healthy): JDM 0.76, against CANDLE 3.15, SAVI 3.13, UAID 2.49, AGS 1.60; healthy groups -0.15 and -0.33. Removing the 28 interferon genes removes about 26 genes from each shared list and leaves most of it (CANDLE 168 to 142, SAVI 154 to 129, UAID 130 to 104). The correlation barely moves because 28 genes are about 0.25% of the genes compared, so that column is a weak test.
- **08, SLE.** SLE against healthy (157 children with SLE, 46 healthy, 13,261 genes, 5,335 with adjusted p < 0.05) matches Anne's run on the same data (157, 46, 13,259, 5,335). Rank correlation with SLE: SAVI 0.46, UAID 0.44, JDM 0.40, CANDLE 0.37, AGS -0.23. These four positive values are close together and should not be ranked. SLE interferon score 2.75 against -0.29 for its healthy children (24 of 28 genes measured).
- **AGS.** `AGS_vs_HC` has 3 batches holding both groups (16 AGS, 62 healthy), so it is not a facility artefact. Its signal is weaker (88 genes up, 114 down at 2-fold) and it is negatively correlated with JDM and SLE. Unexplained; a cell-mix shift is one untested guess. Check: the `gene_name` column of `top_genes_AGS_vs_HC.csv`.
- **Small groups.** `AGS-like` (2 to 3 samples) and `CANDLE-like` (3) appear as their own diagnoses. Notebook 08 skips them by default.

- **09, WGCNA.** RGM network: 113 samples, 27,768 genes, power 12 (the WGCNA FAQ fallback: the scale-free fit never passed R-squared 0.75), 34 modules. JDM-lab network: 112 samples, 12,005 genes, power 12, 29 modules. **Interferon module** in both labs holds all 28 of the 28 interferon genes and none of the NF-kB control genes (RGM `grey60`, 386 genes, median kME 0.88; JDM `tan`, 256 genes, kME 0.91); the two modules share 183 genes (adjusted p about 2e-274), and the RGM one is strongly preserved in the JDM lab. Most other disease-linked RGM modules look like blood cell composition (hubs suggest neutrophil, monocyte, B-cell, eosinophil, NK/T) with modest links (|r| 0.2 to 0.5). Modules with a strong batch link (`steelblue` 0.90, `red` 0.66, `darkorange` 0.63, `orange` 0.50, `brown` 0.42) have mostly non-coding or pseudogene hubs and could not be tested or were not reproduced in the JDM lab, which supports calling them technical. The `batch?` flag in the summary is crude (15 of 34 flagged, many borderline); `module_verdicts_RGM.csv` uses a threshold of 0.4. 27 of 29 JDM modules are preserved in RGM; 8 of 34 RGM modules are strongly preserved, 13 weakly, 7 not, 6 untestable.

**Cautions that apply to all of the above.** Array and RNA-seq are never on one scale; only fold changes against each study's own healthy controls are compared, by rank and direction. A resemblance shared by SLE, JDM, CANDLE, SAVI and UAID may be a general blood signal (interferon, neutrophils) and not specific to any one disease. The corrected table from 04 is for plots and clustering, never for statistical tests.

## Not built yet, and the plan

**WGCNA (notebook 09)** is written, not yet run. It needs the package, installed with `conda install -y -p /home/jovyan/envs/endotypes -c conda-forge r-wgcna` (done), and the R kernel restarted. It follows Anne's notebooks 13 to 15: signed network, power from the scale-free rule (FAQ value if none reaches R-squared 0.90), `minModuleSize` 30, `mergeCutHeight` 0.25, `maxBlockSize` 5000, count filter and PFlog1pPF for the RGM counts. One network per study (RGM from the notebook-03 ComBat-seq counts; the JDM lab from log2 RPKM), then module preservation in both directions and a gene-overlap table that says which RGM module matches which JDM module. Her array modules (GSE65391) are not projected into your data yet.

**Optional side checks.** The broader interferon-gene test (remove the genes most correlated with the interferon score, then redo 07 and 08); the AGS gene-list check above; GSE69529 (43 healthy children) and GSE135779 (PBMC) as further SLE references.

## Folder layout

```
README.md, OUTPUTS.md        this file; what every notebook writes
ipynb/
  nanostring/                01_rcc_parse_qc, 02_append_jdm
  rnaseq/
    nih/                     01 to 10: the run order above (your RGM + JDM-lab data, then SLE and WGCNA)
    public/                  public datasets, run on their own; names are <GSE>_<step>_<what>
      GSE65391_01_download, _02_metadata, _03_normalize         (SLE array; needed by nih/08)
      GSE232381_01_download, _02_metadata, _03_normalise        (public bulk RNA-seq)
      GSE135779_01_download, _02_metadata, _03_pseudobulk, _04_normalise   (SLE single-cell; Anne's version, not needed by nih/10)
src/paths.R                  where files live on Lifebit (everything relative to this folder)
genes/                       the six gene lists (28 interferon response genes and subsets, NF-kB controls, 6-gene score)
env/                         endotypes.yml + setup_endotypes_env.sh (used by step 2); environment.yml (Anne's all-in-one file, for reference)
results/, data/              created on first run: outputs and public downloads (not code; do not commit)
```

## Old names

The RNA-seq chain was renumbered so it runs 01 to 10 with no gap and no "b" notebook. Results text above and `OUTPUTS.md` use the new numbers.

| new (`ipynb/rnaseq/nih/`) | old |
|---|---|
| `01_facility_csvs` | `02_NIH_bulk_RNA_seq_make_facility_csvs` |
| `02_parse_metadata` | `03_NIH_bulk_RNA_seq_parse_metadata` |
| `03_batch_combatseq` | `04_NIH_bulk_RNA_seq_batch_PCA_ComBatSeq` |
| `04_batch_all_sources` | `04b_NIH_JDM_batch_correction_all_sources` |
| `05_differential_expression` … `07_jdm_comparison` | `05_`/`06_`/`07_NIH_bulk_RNA_seq_…` (same numbers) |
| `08_sle_gse65391_vs_jdm`, `09_wgcna_modules`, `10_sle_gse135779_vs_jdm` | same numbers, longer names |

Public: `001-003` are now `GSE65391_01-03`, `004-006` are `GSE232381_01-03`, `007-010` are `GSE135779_01-04`. NanoString: `01_NanoString_RCC_parse_qc` is `01_rcc_parse_qc`, `02_NanoString_append_JDM` is `02_append_jdm`.

Outputs are written to `results/` and public downloads to `data/`, both inside this folder on Lifebit.
