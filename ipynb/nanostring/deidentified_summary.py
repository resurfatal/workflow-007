"""
deidentified_summary.py  -  prints a copy-pasteable, de-identified description of your NanoString files.

Run it on Lifebit (in a notebook cell:  %run /home/jovyan/session_data/workflow-007/ipynb/nanostring/deidentified_summary.py
or paste the whole file into one cell). It READS your files and PRINTS text. It writes nothing.

What it never prints: patient IDs, initials, names, exact dates, or any individual row.
What it prints: column names, how many rows are filled, ID *patterns* (every digit shown as 9, every letter as A),
year ranges of dates, group counts (any group smaller than MIN_CELL is shown as "<5"), and count ranges.
Read the output once before pasting, as you would anything else.
"""
import os, re, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")

# ---- settings ---------------------------------------------------------------
OUT        = "/mnt/filesystems/nanostring_out"
SEARCH     = ["/mnt/filesystems", "/mnt/session_data", "/home/jovyan"]
MASTER_CSV = os.path.join(OUT, "nanostring_master.csv")
WIDE_ALL   = os.path.join(OUT, "nanostring_counts_wide_all_probes.csv")
XLSX_REGEX = r"^JDM IFN scores OADR.*\.xlsx$"
DB_REGEX   = r"^OADR-I-26-007.*clinical.*JDM.*\.csv$"
XLSX = None; DB = None
MIN_CELL = 5                 # groups smaller than this are printed as "<5"
MAX_COLS_LISTED = 80         # long column lists are cut after this many
GENES = ["CXCL10","DDX60","EPSTI1","GBP1","HERC5","HERC6","IFI27","IFI44","IFI44L","IFI6","IFIT1","IFIT2","IFIT3",
         "IFIT5","ISG15","LAMP3","LY6E","MX1","OAS1","OAS2","OAS3","OASL","RSAD2","RTP4","SIGLEC1","SOCS1","SPATS2L","USP18"]
HK = ["ALAS1", "HPRT1", "TBP", "TUBB"]
# columns whose VALUES must never be shown even as group counts (free text / identifiers)
NEVER_SHOW = r"id|initial|name|comment|note|file|mrn|dob|birth|raw|label|cartridge|lane|operator|position|reason"

def pat(x):                  # 015-01-01 -> 999-99-99 ; O-6033 -> A-9999
    return re.sub(r"[A-Za-z]", "A", re.sub(r"\d", "9", str(x).strip()))
def pat_summary(s, k=6):
    s = s.dropna().astype(str); s = s[s.str.strip() != ""]
    if s.empty: return "(empty)"
    vc = s.map(pat).value_counts()
    return "; ".join(f"{p} x{('<%d' % MIN_CELL) if n < MIN_CELL else n}" for p, n in vc.head(k).items())
def cell(n): return f"<{MIN_CELL}" if 0 < n < MIN_CELL else str(n)
def counts(s, k=15):
    vc = s.fillna("(blank)").astype(str).value_counts()
    out = [f"{i}: {cell(n)}" for i, n in vc.head(k).items()]
    if len(vc) > k: out.append(f"... {len(vc)-k} more groups")
    return "; ".join(out)
def to_date(s):
    s = pd.Series(s); num = pd.to_numeric(s, errors="coerce"); d = pd.to_datetime(s, errors="coerce")
    ser = num.between(20000, 80000); d[ser] = pd.to_datetime(num[ser], unit="D", origin="1899-12-30"); return d
def years(s):
    d = to_date(s).dropna()
    return "none parse" if d.empty else f"{int(d.dt.year.min())} to {int(d.dt.year.max())}  (n parsed {len(d)}, distinct years {d.dt.year.nunique()})"
def find_one(pat_):
    hits = []
    for d in SEARCH:
        if os.path.isdir(d):
            for root, _, files in os.walk(d):
                hits += [os.path.join(root, f) for f in files if re.search(pat_, f, re.I) and not f.startswith("~$")]
    return sorted(hits)[0] if hits else None
def head(t): print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)
def collist(df, show_dtype=True):
    rows = []
    for c in df.columns[:MAX_COLS_LISTED]:
        rows.append(f"  {c!r:45} filled {int(df[c].notna().sum()):>5}" + (f"  {df[c].dtype}" if show_dtype else ""))
    if len(df.columns) > MAX_COLS_LISTED: rows.append(f"  ... {len(df.columns)-MAX_COLS_LISTED} more columns not listed")
    print("\n".join(rows))
def count_range(df, cols):
    x = df[cols].apply(pd.to_numeric, errors="coerce")
    q = x.stack().quantile([0, .5, 1]).round(1).tolist()
    whole = bool(np.all(np.isclose(x.stack().dropna() % 1, 0)))
    print(f"  counts across these genes: min {q[0]}, median {q[1]}, max {q[2]}; all whole numbers: {whole}; missing cells: {int(x.isna().sum().sum())}")
    per = x.median().round(0)
    print("  median per gene:", ", ".join(f"{g} {int(v) if pd.notna(v) else 'NA'}" for g, v in per.items()))

# ---- A. RGM master table (notebook 01) ---------------------------------------
head("A. RGM master table (nanostring_master.csv from notebook 01)")
if os.path.exists(MASTER_CSV):
    m = pd.read_csv(MASTER_CSV, low_memory=False)
    print(f"rows {len(m)}, columns {m.shape[1]}")
    print("\ncolumns (name, how many rows filled, type):"); collist(m)
    for c in ["record_type", "diagnosis1", "specimen_type", "qc_pass", "is_non_blood", "needs_review", "study_status", "gene_rlf", "n_probes", "n_housekeeping"]:
        if c in m: print(f"\n{c}: {counts(m[c])}")
    for c in ["sample_date", "run_date"]:
        if c in m: print(f"\n{c} years: {years(m[c])}")
    for c in ["patient_id", "sample_initials", "cartridge", "raw_filename"]:
        if c in m: print(f"\n{c} pattern (digits=9, letters=A): {pat_summary(m[c])}")
    if "patient_id" in m:
        pid = m["patient_id"].replace("", np.nan).dropna()
        print(f"\npatients: {pid.nunique()} distinct; files per patient: " + counts(pid.value_counts().clip(upper=5).astype(int).astype(str) + (" files"), 6))
    if "gene_rlf" in m and "record_type" in m:
        print("\ncodeset (gene_rlf) by record_type  [do RGM files mix codesets? the SOP batch-calibrates Raphaela II vs III]:")
        ct = pd.crosstab(m["gene_rlf"].fillna("(blank)"), m["record_type"])
        print(ct.map(lambda n: cell(int(n))).to_string())
    if "cartridge" in m and "record_type" in m:
        cp = m.groupby("cartridge")["record_type"].agg(lambda s: ", ".join(sorted(set(s.astype(str)))))
        n_c = len(cp); both = int(cp.str.contains("healthy_control") .sum()); stdc = int(cp.str.contains("standard").sum())
        print(f"\ncartridges: {n_c}; cartridges holding a healthy control: {cell(both)}; holding a standard: {cell(stdc)}")
        if "diagnosis1" in m:
            dx_c = m[m["record_type"] == "patient_sample"].groupby("cartridge")["diagnosis1"].nunique()
            print(f"patient-sample cartridges holding >1 diagnosis: {int((dx_c > 1).sum())} of {len(dx_c)}")
    genes_here = [g for g in GENES + HK if g in m]
    print(f"\npanel + housekeeping genes present as columns: {len(genes_here)} of {len(GENES)+len(HK)}; missing: {[g for g in GENES+HK if g not in m]}")
    if genes_here: count_range(m, genes_here)
    qc = [c for c in ["binding_density", "fov_ratio", "pos_control_geomean", "neg_control_mean", "limit_of_detection", "hk_geomean", "total_counts"] if c in m]
    if qc: print("\nQC columns (min / median / max):\n" + "\n".join(f"  {c}: " + " / ".join(f"{v:.2f}" for v in pd.to_numeric(m[c], errors='coerce').quantile([0,.5,1])) for c in qc))
else:
    print("not found:", MASTER_CSV)

head("A2. Does the all-probes file contain the control probes? (needed for positive-control scaling)")
if os.path.exists(WIDE_ALL):
    w = pd.read_csv(WIDE_ALL, nrows=3)
    ctrl = [c for c in w.columns if re.match(r"(POS|NEG)_", str(c), re.I)]
    print(f"{w.shape[1]} columns; control probe columns: {ctrl if ctrl else 'NONE'}")
    others = [c for c in w.columns if c not in GENES + HK and not re.match(r"(POS|NEG)_", str(c), re.I)]
    print("other non-gene columns (first 25 names):", others[:25])
else:
    print("not found:", WIDE_ALL, "(it is only written when INCLUDE_ALL_PROBES is on in notebook 01; the counts_long file would also do)")
    lg = os.path.join(OUT, "nanostring_counts_long.csv.gz")
    if os.path.exists(lg):
        l = pd.read_csv(lg, nrows=200000)
        print("long file exists; columns:", list(l.columns)[:12])
        cc = [c for c in l.columns if re.search("class", c, re.I)]
        if cc: print("code classes:", counts(l[cc[0]]))

# ---- B. JDM workbook ---------------------------------------------------------
head("B. JDM workbook")
XLSX = XLSX or find_one(XLSX_REGEX)
if XLSX:
    xl = pd.ExcelFile(XLSX); print("sheets:", xl.sheet_names)
    for sh in xl.sheet_names:
        raw = pd.read_excel(XLSX, sheet_name=sh, header=None, dtype=str).dropna(how="all")
        if raw.empty: print(f"\n--- sheet {sh!r}: empty"); continue
        hd = raw.head(15)
        hits = [sum(str(x).strip() in set(GENES + HK) for x in hd.iloc[i]) for i in range(len(hd))]
        h = int(np.argmax(hits)); hdr = [str(x).strip() if pd.notna(x) else f"unnamed_{j}" for j, x in enumerate(hd.iloc[h])]
        df = raw.iloc[h + 1:].copy(); df.columns = pd.Index(hdr).where(~pd.Index(hdr).duplicated(), pd.Index([f"{x}__dup{i}" for i,x in enumerate(hdr)]))
        print(f"\n--- sheet {sh!r}: {len(df)} data rows, header at row {h} (genes found in header: {hits[h]} of {len(GENES)+len(HK)})")
        nongene = [c for c in df.columns if c not in GENES + HK]
        print("non-gene columns:", nongene[:40])
        for c in nongene:
            if re.search(r"sample|patient|initial|study", c, re.I): print(f"  pattern of {c!r}: {pat_summary(df[c])}")
            if re.search(r"date", c, re.I): print(f"  years in {c!r}: {years(df[c])}")
        g_here = [g for g in GENES + HK if g in df.columns]
        if g_here: count_range(df, g_here)
        score_cols = [c for c in nongene if re.search(r"geomean|std|score|ifn", c, re.I)]
        if score_cols:
            print("  workbook's own score columns (min / median / max):")
            for c in score_cols:
                v = pd.to_numeric(df[c], errors="coerce")
                if v.notna().sum(): print(f"    {c}: n {int(v.notna().sum())}, " + " / ".join(f"{x:.1f}" for x in v.quantile([0, .5, 1])))
        sid = [c for c in df.columns if re.search("study sample", c, re.I)]
        if sid:
            pt = df[sid[0]].astype(str).str.strip().str.split("-").str[0]
            if "JDM_PATIENTS" not in globals(): JDM_PATIENTS = set(pt.str.upper())
            vc = pt.value_counts()
            print(f"  patients (text before first '-'): {pt.nunique()}; samples per patient: " + counts(vc.clip(upper=5).astype(str), 6))
else:
    print("workbook not found")

# ---- C. clinical database ----------------------------------------------------
head("C. Clinical database (column names and fill rates only; no values)")
DB = DB or find_one(DB_REGEX)
if DB:
    d = pd.read_csv(DB, dtype=str, low_memory=False)
    print(f"rows {len(d)}, columns {d.shape[1]}")
    key = [c for c in d.columns if re.search(r"date|visit|nano|activity|ACT|DAS|MMT|CK|MSA|P155|MDA5|MJ|DX|duration|RNA", c, re.I)]
    print(f"\ncolumns that matter here ({len(key)} of {d.shape[1]}):")
    collist(d[key], show_dtype=False)
    print("\nall other column names (no fill rates):", [c for c in d.columns if c not in key][:150])
    for c in key:
        if re.search(r"nano|RNA", c, re.I): print(f"\npattern of {c!r}: {pat_summary(d[c])}")
        if re.search(r"date", c, re.I) and not re.search(r"birth|dob", c, re.I): print(f"years in {c!r}: {years(d[c])}")
    for c in ["DX_PROTOCOL", "Basic_MSA", "P155", "MDA5", "MJ"]:
        if c in d: print(f"\n{c}: {counts(d[c])}")
else:
    print("clinical database not found")

# ---- D. overlap, counts only -------------------------------------------------
head("D. Overlap (counts only)")
try:
    if XLSX and DB and "JDM_PATIENTS" in globals():
        norm = lambda s: s.astype(str).str.strip().str.split("-").str[0].str.upper()
        jp = JDM_PATIENTS
        nc = [c for c in d.columns if c.lower() == "nano id"]
        if nc:
            dp = set(norm(d[nc[0]].dropna()))
            print(f"JDM patients {len(jp)}; in clinical DB {len(jp & dp)}; not in DB {len(jp - dp)}")
except Exception as e:
    print("overlap skipped:", e)
print("\nDone. Copy everything above this line.")
