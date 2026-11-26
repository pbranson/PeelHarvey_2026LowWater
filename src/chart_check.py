"""Datum check against DoT nautical charts WA848 / WA913 / WA001 tidal-plane tables.

The chart tables give MSL above the chart's sounding/tidal datum (CD) and tie that CD
to a benchmark. Where the gauge header ties its own datum to the *same* benchmark we
can convert chart MSL to the gauge datum and compare it with the gauge's observed
mean level. The chart tidal-plane epoch is not stated, so gauge means are shown for
several windows; sea-level rise (~3-6 mm/yr here) changes MSL by ~3-6 cm/decade.

Run after run_analysis.py (uses data/processed/hourly_raw.pkl).
"""
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import datums  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# Chart tidal planes (m above chart datum), read from the chart PDFs
CHART = {
    "Mandurah Ocean Entrance (WA848)": dict(HAT=1.26, MHHW=0.91, MLHW=0.79, MSL=0.77, MHLW=0.70, MLLW=0.62, LAT=0.35,
                                             datum="2.060 m below Tidal BM 1971"),
    "Peel Inlet (WA848)": dict(HAT=0.93, MHHW=0.72, MLHW=0.56, MSL=0.56, MHLW=0.48, MLLW=0.43, LAT=0.31,
                               datum="1.216 m below BM A895"),
    "Dawesville Ocean Entrance (WA848)": dict(HAT=1.24, MHHW=0.81, MLHW=0.62, MSL=0.62, MHLW=0.52, MLLW=0.42, LAT=0.17,
                                              datum="3.491 m below BM DMH 102"),
    "Harvey Estuary (WA848)": dict(HAT=1.02, MHHW=0.71, MLHW=0.58, MSL=0.59, MHLW=0.50, MLLW=0.45, LAT=0.30,
                                   datum="1.522 m below SSM PHS 6"),
    "Ocean (WA913)": dict(HAT=1.32, MHHW=0.96, MLHW=0.70, MSL=0.73, MHLW=0.72, MLLW=0.52, LAT=0.14, datum="n/a"),
    "Peel Inlet (WA913)": dict(HAT=0.93, MHHW=0.72, MLHW=0.56, MSL=0.56, MHLW=0.48, MLLW=0.43, LAT=0.31, datum="n/a"),
    "Mandurah Estuary (WA913)": dict(HAT=1.01, MHHW=0.74, MLHW=0.60, MSL=0.60, MHLW=0.52, MLLW=0.45, LAT=0.27, datum="n/a"),
    "Fremantle (WA001)": dict(HAT=1.43, MHHW=1.18, MLHW=1.08, MSL=0.84, MHLW=0.61, MLLW=0.50, LAT=0.29,
                              datum="3.088 m below BM 53A"),
}

# gauge datum vs chart datum where they share a benchmark: chart CD = gauge datum + delta
SHARED_BM = {
    # Peel gauge (2017+ segment / 2026 file): CD 1.201 m below A895; chart 1.216 m below A895
    "PLPEE01": dict(chart="Peel Inlet (WA848)", delta=-(1.216 - 1.201)),
    # Mandurah Fishermans Jetty: CD 2.065 m below TBM 1971; chart 2.060 m below TBM 1971
    "MHMAN01": dict(chart="Mandurah Ocean Entrance (WA848)", delta=+(2.065 - 2.060)),
}
WINDOWS = [("1991", "1999"), ("1996", "2005"), ("2001", "2010"), ("2007", "2016"), ("2017", "2025")]


def main():
    raw = pd.read_pickle(os.path.join(ROOT, "data", "processed", "hourly_raw.pkl")).drop(columns="HILL_P")
    M0 = datums.monthly_means(raw)
    steps = datums.estimate_steps(M0)
    hom = datums.apply_steps(raw, steps)
    M = datums.monthly_means(hom)
    off = datums.ahd_offsets(M)
    A = M.copy()
    for st, o in off.items():
        A[st] = A[st] + o
    ocean_ahd = A["HILLARYS"]

    rows = []
    for st, info in SHARED_BM.items():
        ch = CHART[info["chart"]]
        chart_msl_gauge = ch["MSL"] + info["delta"]  # chart MSL expressed on gauge datum
        for a, b in WINDOWS:
            g = M[st].loc[a:b]
            if g.notna().sum() < 36:
                continue
            rows.append(dict(gauge=st, chart=info["chart"], window=f"{a}-{b}",
                             gauge_mean_m=round(g.mean(), 3), chart_MSL_on_gauge_datum_m=round(chart_msl_gauge, 3),
                             gauge_minus_chart_cm=round(100 * (g.mean() - chart_msl_gauge), 1)))
    shared = pd.DataFrame(rows)

    # Implied MSL in AHD at the chart sites, using this analysis's AHD reductions,
    # compared with ocean MSL (Hillarys AHD) in the same window
    imp = []
    for st, chart in [("PLPEE01", "Peel Inlet (WA848)"), ("MHMAN01", "Mandurah Ocean Entrance (WA848)"),
                      ("DVHAR01", "Harvey Estuary (WA848)")]:
        for a, b in WINDOWS:
            g = A[st].loc[a:b]
            if g.notna().sum() < 36:
                continue
            imp.append(dict(gauge=st, window=f"{a}-{b}", gauge_mean_mAHD=round(g.mean(), 3),
                            hillarys_mean_mAHD=round(ocean_ahd.loc[a:b].mean(), 3),
                            estuary_minus_ocean_cm=round(100 * (g.mean() - ocean_ahd.loc[a:b].mean()), 1)))
    implied = pd.DataFrame(imp)

    # Harvey & Caddadup have no shared benchmark: infer their gauge datum vs chart CD
    harv = []
    for a, b in WINDOWS:
        g = M["DVHAR01"].loc[a:b]
        if g.notna().sum() >= 36:
            harv.append(dict(window=f"{a}-{b}", gauge_mean_m=round(g.mean(), 3),
                             implied_chartCD_minus_gaugeDatum_cm=round(100 * (g.mean() - CHART["Harvey Estuary (WA848)"]["MSL"]), 1)))

    # Chart MSL differences between sites vs gauge-observed mean differences (datum-free check)
    # e.g. Peel - Mandurah estuary mean level difference cannot be checked without common datum: report chart values
    out = dict(chart_tidal_planes=CHART, shared_benchmark_check=shared.to_dict("records"),
               implied_msl_ahd=implied.to_dict("records"), harvey_gauge_vs_chart=harv)
    with open(os.path.join(ROOT, "outputs", "chart_datum_check.json"), "w") as f:
        json.dump(out, f, indent=2)
    pd.set_option("display.width", 200)
    print(shared.to_string(index=False)); print()
    print(implied.to_string(index=False)); print()
    print(pd.DataFrame(harv).to_string(index=False))
    # raw (un-homogenised) Peel segment means vs chart, to show which segments share the chart CD
    print("\nPeel raw segment means vs chart MSL on gauge datum (0.545 m):")
    r = datums.monthly_means(raw)["PLPEE01"]
    for a, b in [("1987", "1994"), ("1995", "1999"), ("2000", "2016"), ("2017", "2025")]:
        print(f"  {a}-{b}: {r.loc[a:b].mean():.3f}")
    r = datums.monthly_means(raw)["MHMAN01"]
    print("Mandurah FJ raw segment means vs chart MSL on gauge datum (0.775 m):")
    for a, b in [("1991", "2000"), ("2002", "2006")]:
        print(f"  {a}-{b}: {r.loc[a:b].mean():.3f}")


if __name__ == "__main__":
    main()
