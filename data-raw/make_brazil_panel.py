"""Build the four-column Brazil property-tax panel bundled with PanelAdequacy.

Usage:
    python make_brazil_panel.py --source Christensen_Garfias_2021_JOP.rds \
        --output brazil_property_tax_panel.csv

The transformation matches Paper C: retain complete 2004--2015 municipal
paths with nonmissing outcome and treatment, then exclude municipalities
already treated in 2004. The output uses consecutive unit/period codes and a
blank first-treatment value for never-treated municipalities.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import pyreadr


def read_frame(path: Path) -> pd.DataFrame:
    objects = [value for value in pyreadr.read_r(str(path)).values()
               if isinstance(value, pd.DataFrame)]
    if len(objects) != 1:
        raise ValueError(f"expected one data frame in {path}, found {len(objects)}")
    return objects[0]


def build(source: Path) -> tuple[pd.DataFrame, dict[str, int]]:
    frame = read_frame(source)
    required = ["c6_ibge", "year", "cad_update", "logiptu"]
    available = frame.dropna(subset=required).copy()
    years = sorted(frame["year"].dropna().unique())
    periods = len(years)
    counts = available.groupby("c6_ibge")["year"].nunique()
    complete_ids = counts[counts == periods].index
    complete = available[available["c6_ibge"].isin(complete_ids)].copy()
    complete.sort_values(["c6_ibge", "year"], inplace=True)
    initially_treated = complete.groupby("c6_ibge")["cad_update"].first() != 0
    eligible_ids = initially_treated[~initially_treated].index
    panel = complete[complete["c6_ibge"].isin(eligible_ids)].copy()

    unit_map = {value: i + 1 for i, value in enumerate(sorted(eligible_ids))}
    time_map = {value: i + 1 for i, value in enumerate(years)}
    panel["uid"] = panel["c6_ibge"].map(unit_map).astype(int)
    panel["tid"] = panel["year"].map(time_map).astype(int)

    treated = panel[panel["cad_update"] != 0]
    first = treated.groupby("uid")["tid"].min()
    panel["ft"] = panel["uid"].map(first).astype("Int64")
    output = panel[["uid", "tid", "ft", "logiptu"]].rename(
        columns={"logiptu": "y"})

    if len(output) != 34_080 or output["uid"].nunique() != 2_840:
        raise ValueError("Brazil extract does not match the 2,840 x 12 audit sample")
    treatment = (output["ft"].notna() &
                 (output["tid"] >= output["ft"])).astype(int)
    if (treatment.groupby(output["uid"]).diff().fillna(0) < 0).any():
        raise ValueError("treatment is not absorbing")
    audit = {
        "raw_rows": len(frame),
        "raw_units": frame["c6_ibge"].nunique(),
        "periods": periods,
        "complete_units": len(complete_ids),
        "initially_treated_complete_units": int(initially_treated.sum()),
        "eligible_units": len(eligible_ids),
        "eligible_rows": len(output),
    }
    return output, audit


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output, audit = build(args.source)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False, na_rep="")
    print(audit)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
