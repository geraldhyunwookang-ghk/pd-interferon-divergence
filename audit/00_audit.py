"""
00_audit.py — pre-freeze donor audit (plan v2.0, section 4 and Appendix A).
Reads ONLY donor metadata and nucleus-annotation files. No expression values.
Run from a folder containing the metadata files listed in audit/sources.csv.
Writes audit/donor_table.csv (one row per donor per cohort).
"""
import gzip, os, tarfile, openpyxl, pandas as pd

def first(*names):
    for n in names:
        if os.path.exists(n): return n
    raise FileNotFoundError(names)

def series_matrix(path):
    rows = {}
    for line in gzip.open(path, "rt"):
        if line.startswith("!Sample_title") or line.startswith("!Sample_characteristics"):
            f = [x.strip('"') for x in line.rstrip("\n").split("\t")][1:]
            k = "title" if line.startswith("!Sample_title") else f[0].split(":")[0]
            rows.setdefault(k, [v.split(": ", 1)[-1] for v in f])
    return pd.DataFrame(rows)

out = []

# ---- Kamath (SCP1768) ------------------------------------------------------
meta = pd.read_csv(gzip.open(first("METADATA_PD.tsv.gz", "METADATA_PD.tsv"), "rt"), sep="\t", skiprows=[1], dtype=str,
                   usecols=["NAME", "donor_id", "Status", "sex", "Donor_Age", "Donor_PMI"])
key = meta.set_index("NAME")
donors = meta.groupby("donor_id").first()
def count(f, filt=None):
    d = pd.read_csv(f, sep="\t", skiprows=[1])
    if filt is not None: d = d[filt(d)]
    return d.join(key, on="NAME").groupby("donor_id").size()
da = pd.read_csv("da_UMAP.tsv", sep="\t", skiprows=[1]).join(key, on="NAME")
k = pd.DataFrame({
    "n_DA": da.groupby("donor_id").size(),
    "n_DA_SOX6": da[da.Cell_Type.str.startswith("SOX6")].groupby("donor_id").size(),
    "n_DA_CALB1": da[da.Cell_Type.str.startswith("CALB1")].groupby("donor_id").size(),
    "n_microglia": count("mg_UMAP.tsv", lambda d: ~d.Cell_Type.str.startswith("Macro")),
    "n_microglia_incl_macro": count("mg_UMAP.tsv"),
    "n_astro": count("astro_UMAP.tsv"), "n_olig": count("olig_UMAP.tsv"),
    "n_opc": count("opc_UMAP.tsv"), "n_nonDA": count("nonda_UMAP.tsv"),
}).fillna(0).astype(int)
k = k.join(donors[["Status", "sex", "Donor_Age", "Donor_PMI"]])
k = k.rename(columns={"Status": "diagnosis", "Donor_Age": "age", "Donor_PMI": "pmi_h"})
k["diagnosis"] = k.diagnosis.replace({"Ctrl": "Control"})
k["cohort"] = "Kamath"; k["donor"] = k.index
out.append(k.reset_index(drop=True))

# ---- Smajic (GSE157783) ----------------------------------------------------
sm = series_matrix("GSE157783_series_matrix.txt.gz")
if not os.path.exists("IPDCO_hg_midbrain_cell.tsv"):
    tarfile.open(first("GSE157783_IPDCO_hg_midbrain_cell.tar.gz")).extractall()
cell = pd.read_csv("IPDCO_hg_midbrain_cell.tsv", sep="\t")
sm["donor"] = sm.title.str.replace("^IPD", "PD", regex=True)   # PDk <-> IPDk (plan 4.3)
nmg = cell[cell.cell_ontology == "Microglia"].groupby("patient").size()
nda = cell[cell.cell_ontology == "DaNs"].groupby("patient").size()
s = pd.DataFrame({"cohort": "Smajic", "donor": sm.donor, "geo_title": sm.title,
    "diagnosis": ["PD" if "Parkinson" in x else "Control" for x in sm["disease state"]],
    "age": sm["age at death"], "sex": sm.Sex, "pmi_h": sm["pmi (h)"]})
s["n_microglia"] = s.donor.map(nmg).fillna(0).astype(int)
s["n_DA"] = s.donor.map(nda).fillna(0).astype(int)
out.append(s)

# ---- Martirosyan (GSE243639) -----------------------------------------------
mm = series_matrix("GSE243639_series_matrix.txt.gz")
mm["donor"] = mm.title.str.extract(r"(s_\d+)")[0].str.replace("_", ".")
ws = openpyxl.load_workbook("13024_2023_699_MOESM2_ESM.xlsx", data_only=True)["Nuclei counts per donor"]
st = pd.DataFrame([(str(r[2]), int(r[7] or 0), int(r[13] or 0))
                   for r in ws.iter_rows(min_row=5, max_row=33, values_only=True) if r[0] is not None],
                  columns=["donor", "n_microglia", "n_DA"]).set_index("donor")
m = pd.DataFrame({"cohort": "Martirosyan", "donor": mm.donor,
    "diagnosis": ["PD" if "Parkinson" in x else "Control" for x in mm["clinical diagnosis"]],
    "age": mm.age, "sex": mm.Sex, "pmi_h": mm["pmi hours"], "rin": mm["rin measure"]})
m = m.join(st, on="donor")
out.append(m)

cols = ["cohort", "donor", "geo_title", "diagnosis", "age", "sex", "pmi_h", "rin",
        "n_DA", "n_DA_SOX6", "n_DA_CALB1", "n_microglia", "n_microglia_incl_macro",
        "n_astro", "n_olig", "n_opc", "n_nonDA"]
t = pd.concat(out, ignore_index=True).reindex(columns=cols)
t["sex"] = t.sex.str.lower().replace({"f": "female", "m": "male"})
t.to_csv("audit/donor_table.csv", index=False)
print(t.groupby(["cohort", "diagnosis"]).size().to_string())
