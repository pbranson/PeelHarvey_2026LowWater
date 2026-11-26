"""Datum homogenisation and AHD reduction, estimated on monthly means.

Reference: BoM ABSLMP Hillarys (IDO71012), AHD = gauge zero + 0.763 m (ABSLMP
station table). Fremantle DoT FFFBH01 is on LWM Fremantle 1949 = AHD - 0.756 m;
the two agree to within ~2 cm (see README), which confirms the reduction.

Datum *steps* (instrument/datum changes) are removed by comparing the
station-minus-Hillarys monthly difference in windows either side of a known
break. Each record is aligned to its most recent segment. Slow drift common to
all Peel stations relative to Hillarys (~1.5-2 mm/yr, consistent with
documented Perth Basin subsidence at Hillarys) is NOT removed.

Stations without a published AHD relationship are reduced to AHD by matching
monthly means to the reference in a common period:
  * open-coast gauges (Cape Bouvard, Mandurah Ocean Marina): all months;
  * estuary gauges (Peel Inlet, Caddadup, Harvey): low-runoff months (Jan-Mar)
    against Mandurah, assuming negligible estuary-ocean mean-level gradient in
    late summer. This is an *approximation*; replace with surveyed benchmark
    AHD heights (A895, MAN TIDAL 02, Florida 2) when available.
"""
import pandas as pd

HILLARYS_AHD = 0.763   # m, gauge zero -> AHD (ABSLMP)
FREO_LWMF_AHD = 0.756  # m, LWM Fremantle 1949 is 0.756 m below AHD

# (station, break date, months to exclude after the break, window months, months-of-year to use)
# Peel/Harvey breaks use Nov-Apr only so the post-Dawesville-Cut (Apr 1994) change
# in winter setup is not mistaken for a datum step.
BREAKS = [
    ("PLPEE01", "1995-01-01", 0, 12, [11, 12, 1, 2, 3, 4]),
    ("DVHAR01", "1997-01-01", 0, 12, [11, 12, 1, 2, 3, 4]),
    ("DVBVD01", "2013-02-01", 0, 24, None),   # step across 2013-2015 outage
    ("FFFBH01", "2013-01-01", 0, 24, None),
]
TRIM_START = {"DVBVD01": "1995-01-01"}  # early sensor settling 1993-94


def monthly_means(hourly: pd.DataFrame, min_days: int = 20) -> pd.DataFrame:
    cnt = hourly.resample("MS").count()
    return hourly.resample("MS").mean().where(cnt >= 24 * min_days)


def estimate_steps(M: pd.DataFrame, ref: str = "HILLARYS") -> list:
    steps = []
    for st, date, gap, win, months in BREAKS:
        d = (M[st] - M[ref]).dropna()
        t0 = pd.Timestamp(date)
        before = d[(d.index < t0) & (d.index >= t0 - pd.DateOffset(months=win))]
        # for a break inside an outage, the window starts at the first valid month after it
        t1 = d.index[d.index >= t0 + pd.DateOffset(months=gap)][0]
        after = d[(d.index >= t1) & (d.index < t1 + pd.DateOffset(months=win))]
        if months:
            before = before[before.index.month.isin(months)]
            after = after[after.index.month.isin(months)]
        steps.append(dict(station=st, date=date, step=after.mean() - before.mean(),
                          n_before=len(before), n_after=len(after)))
    return steps


def apply_steps(df: pd.DataFrame, steps: list) -> pd.DataFrame:
    """Shift data before each break by the step so earlier data align with later data.
    Works on any time-indexed frame (hourly or monthly)."""
    out = df.copy()
    for s in steps:
        m = out.index < pd.Timestamp(s["date"])
        out.loc[m, s["station"]] += s["step"]
    for st, t in TRIM_START.items():
        if st in out:
            out.loc[out.index < pd.Timestamp(t), st] = float("nan")
    return out


def ahd_offsets(M: pd.DataFrame, period=("2015-01-01", "2025-12-31")) -> dict:
    """Return the constant (m) to ADD to each homogenised station to express it in AHD."""
    p = M.loc[period[0]:period[1]]
    ref = p["HILLARYS"] - HILLARYS_AHD
    off = {"HILLARYS": -HILLARYS_AHD, "FFFBH01": -FREO_LWMF_AHD}
    for st in ("DVBVD01", "MHMAN02"):
        off[st] = (ref - p[st]).mean()
    man = p["MHMAN02"] + off["MHMAN02"]
    summer = p.index.month.isin([1, 2, 3])
    for st in ("PLPEE01", "DVCAD01", "DVHAR01"):
        off[st] = (man - p[st])[summer].mean()
    # Mandurah Fishermans Jetty (CD 2.065 m below TBM 1971): match to Hillarys 2002-2006
    q = M.loc["2002":"2006"]
    off["MHMAN01"] = ((q["HILLARYS"] - HILLARYS_AHD) - q["MHMAN01"]).mean()
    return off
