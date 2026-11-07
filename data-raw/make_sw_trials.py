"""
Build the two stepped-wedge trial extracts bundled with panelcert / PanelAdequacy.

    python data-raw/make_sw_trials.py HAINES_S2.xlsx DOTS99_S1.xlsx OUTDIR

Inputs (public PLOS Medicine Supporting Information, CC BY 4.0):
  Haines et al. (2017), PLOS Med 14(10):e1002412, S2 Data, doi:10.1371/journal.pmed.1002412.s002
  Cattamanchi et al. (2021), PLOS Med 18(5):e1003628, S1 Data, doi:10.1371/journal.pmed.1003628.s002

Outputs:
  haines2017_trial1.csv  Trial 1 analysis extract (study1 == 1): one row per patient episode with
                         ward, site, period (sw_step), start (first treated period of the ward),
                         trt (no_we_exposure: weekend allied-health service withdrawn), acute_los.
                         No other columns (in particular no identifiers or dates) are retained.
  dots99_uganda.csv      the complete S1 Data sheet (all 2,224 rows, 13 columns) with snake_case column
                         names. Values are unchanged except that commas inside text labels are replaced by
                         semicolons (e.g. "P-BC (Bacteriologically confirmed; pulmonary TB)"), so the file
                         can be read by the dependency-free CSV readers of both packages.
"""
import hashlib
import os
import sys

import pandas as pd


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(haines, dots, outdir):
    h = pd.read_excel(haines)
    h = h[h.study1 == 1]
    out = pd.DataFrame({
        "ward": h.index_ward.astype(int),
        "site": h.site.astype(int) if pd.api.types.is_numeric_dtype(h.site) else h.site.astype(str),
        "period": h.sw_step.astype(int),
        "trt": h.no_we_exposure.astype(int),
        "acute_los": h.acute_los.astype(float),
    })
    start = out[out.trt == 1].groupby("ward").period.min()
    out.insert(3, "start", out.ward.map(start).astype(int))
    assert ((out.period >= out.start) == (out.trt == 1)).all()
    out = out.sort_values(["ward", "period"], kind="mergesort").reset_index(drop=True)
    p1 = os.path.join(outdir, "haines2017_trial1.csv")
    out.to_csv(p1, index=False, float_format="%.10g")

    u = pd.read_excel(dots)
    u.columns = ["block", "facility", "facility_type", "month", "period", "sex", "age", "hiv",
                 "disease_class", "retreatment", "treated_successfully", "persisting", "not_ltfu"]
    for col in u.columns:
        if u[col].dtype == object or str(u[col].dtype).startswith("str"):
            u[col] = u[col].astype(str).str.replace(",", ";", regex=False)   # keep the CSV comma-safe
    p2 = os.path.join(outdir, "dots99_uganda.csv")
    u.to_csv(p2, index=False)

    for src, dst in ((haines, p1), (dots, p2)):
        print(f"{os.path.basename(src)}  SHA-256 {sha256(src)}")
        print(f"  -> {os.path.basename(dst)}  SHA-256 {sha256(dst)}  rows {sum(1 for _ in open(dst)) - 1}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
