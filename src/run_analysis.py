"""Peel-Harvey subtidal water level analysis and correlation with the SOI.

Pipeline
  1. Load DoT WA gauges (5/15-min, AWST), BoM ABSLMP Hillarys (hourly, UTC), SOI.
  2. Hourly grid (UTC), short-gap interpolation (<= 6 h).
  3. Datum steps removed and AHD reduction (src/datums.py), on monthly means.
  4. 40-h low-pass Butterworth (4th order, zero-phase filtfilt) -> subtidal level.
  5. Monthly means of subtidal level, anomalies vs 1996-2025 monthly climatology.
  6. Lagged correlation with BoM monthly SOI (raw and 5-month running means),
     significance from effective degrees of freedom (Bretherton et al. 1999).

Run:  python src/run_analysis.py      (writes outputs/ and figures/)
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
import datums  # noqa: E402
import load_data as L  # noqa: E402
from analysis_core import butter_lowpass, load_soi_bom, to_hourly  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "outputs")
FIG = os.path.join(ROOT, "figures")
os.makedirs(OUT, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

NAMES = {
    "HILLARYS": "Hillarys (ocean, BoM)",
    "FFFBH01": "Fremantle (ocean, DoT)",
    "DVBVD01": "Cape Bouvard (ocean)",
    "MHMAN02": "Mandurah Ocean Marina",
    "PLPEE01": "Peel Inlet",
    "DVCAD01": "Caddadup",
    "DVHAR01": "Harvey Estuary",
}
ESTUARY = ["MHMAN02", "PLPEE01", "DVCAD01", "DVHAR01"]
OCEAN = ["HILLARYS", "FFFBH01", "DVBVD01"]
# dataviz default categorical palette, fixed order (validated, light mode)
COL = dict(zip(["HILLARYS", "DVBVD01", "MHMAN02", "PLPEE01", "DVCAD01", "DVHAR01"],
               ["#257de0", "#f47230", "#00b894", "#f2a600", "#ee7fac", "#008a17"]))
COL["FFFBH01"] = "#6a4ce2"
INK, INK2, GRID = "#191919", "#666666", "#e6e6e6"
CLIM = ("1996-01-01", "2025-12-31")   # post Dawesville Cut (opened Apr 1994)

plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
    "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.2,
    "figure.dpi": 150, "savefig.bbox": "tight", "legend.frameon": False,
})


# ----------------------------------------------------------------------------- load
def load_all() -> pd.DataFrame:
    cache = os.path.join(ROOT, "data", "processed", "hourly_raw.pkl")
    if os.path.exists(cache):
        return pd.read_pickle(cache)
    H = {st: to_hourly(L.load_dot(st)) for st in ["FFFBH01", "DVBVD01", "MHMAN02", "MHMAN01",
                                                   "PLPEE01", "DVCAD01", "DVHAR01"]}
    H["HILLARYS"] = to_hourly(L.load_hillarys())
    H["HILL_P"] = to_hourly(L.load_hillarys_pressure())
    df = pd.DataFrame(H)
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    df.to_pickle(cache)
    return df


def eff_dof(x, y):
    """Bretherton et al. (1999) effective sample size for correlated series."""
    x, y = np.asarray(x), np.asarray(y)
    n = len(x)
    r1x = np.corrcoef(x[:-1], x[1:])[0, 1]
    r1y = np.corrcoef(y[:-1], y[1:])[0, 1]
    return max(3.0, n * (1 - r1x * r1y) / (1 + r1x * r1y))


def corr_sig(x, y):
    m = x.notna() & y.notna()
    a, b = x[m], y[m]
    r = np.corrcoef(a, b)[0, 1]
    ne = eff_dof(a, b)
    t = r * np.sqrt((ne - 2) / max(1e-9, 1 - r * r))
    p = 2 * stats.t.sf(abs(t), ne - 2)
    return r, p, ne, int(m.sum())


def main():
    raw = load_all()
    pres = raw.pop("HILL_P")

    # ---- datums (monthly), applied to hourly
    M0 = datums.monthly_means(raw)
    steps = datums.estimate_steps(M0)
    hom = datums.apply_steps(raw, steps)
    Mh = datums.monthly_means(hom)
    off = datums.ahd_offsets(Mh)
    ahd = hom.copy()
    for st, o in off.items():
        ahd[st] = ahd[st] + o
    # Mandurah: splice Fishermans Jetty (1990-2007) ahead of Ocean Marina (2007-)
    ahd["MHMAN02"] = ahd["MHMAN02"].combine_first(ahd["MHMAN01"])
    ahd = ahd.drop(columns="MHMAN01")

    # ---- 40-h Butterworth low-pass
    sub = pd.DataFrame({c: butter_lowpass(ahd[c], cutoff_h=40.0, order=4) for c in ahd})
    # inverse-barometer correction for Hillarys (1 hPa ~ 0.993 cm)
    ib = -(butter_lowpass(pres, 40.0) - 1013.25) * 0.00993
    sub_ib = sub["HILLARYS"] - ib

    sub.to_pickle(os.path.join(ROOT, "data", "processed", "subtidal_hourly.pkl"))
    daily = sub.resample("D").mean().where(sub.resample("D").count() >= 18)
    daily.round(4).to_csv(os.path.join(OUT, "subtidal_daily_mAHD.csv"), index_label="date_utc")

    monthly = sub.resample("MS").mean().where(sub.resample("MS").count() >= 24 * 20)
    monthly["HILLARYS_IB"] = sub_ib.resample("MS").mean().where(sub_ib.resample("MS").count() >= 480)
    clim = monthly.loc[CLIM[0]:CLIM[1]].groupby(monthly.loc[CLIM[0]:CLIM[1]].index.month).mean()
    anom = monthly - clim.reindex(monthly.index.month).set_index(monthly.index)
    # linear trend fitted to 1996-2025 deseasonalised anomalies, removed (extrapolated into 2026)
    trend = {}
    t_yr = (anom.index - pd.Timestamp("2010-01-01")).days / 365.25
    anom_dt = anom.copy()
    for c in anom:
        m = anom[c].notna() & (anom.index >= CLIM[0]) & (anom.index <= CLIM[1])
        if m.sum() < 60:
            continue
        k = np.polyfit(t_yr[m], anom[c][m], 1)
        trend[c] = k[0] * 1000  # mm/yr
        anom_dt[c] = anom[c] - np.polyval(k, t_yr)
    # estuary minus local ocean (Cape Bouvard): isolates estuary-specific (runoff, wind) signal
    for c in ["PLPEE01", "DVCAD01", "DVHAR01", "MHMAN02"]:
        anom_dt[c + "-DVBVD01"] = anom_dt[c] - anom_dt["DVBVD01"]
    monthly.round(4).to_csv(os.path.join(OUT, "subtidal_monthly_mAHD.csv"), index_label="month")
    anom.round(4).to_csv(os.path.join(OUT, "subtidal_monthly_anomaly_m.csv"), index_label="month")
    anom_dt.round(4).to_csv(os.path.join(OUT, "subtidal_monthly_anomaly_detrended_m.csv"), index_label="month")

    soi = load_soi_bom(os.path.join(ROOT, "data", "soi", "soi_bom_monthly.txt"))
    soi = soi.reindex(anom.index)

    # ---- correlations
    rows = []
    lags = range(-3, 13)
    cols = ["HILLARYS", "HILLARYS_IB", "FFFBH01", "DVBVD01", "MHMAN02", "PLPEE01", "DVCAD01", "DVHAR01"]
    for kind, A in (("raw", anom), ("detrended", anom_dt)):
      for smooth in (1, 5):
        a_s = A.rolling(smooth, center=True, min_periods=smooth).mean()
        s_s = soi.rolling(smooth, center=True, min_periods=smooth).mean()
        for st in cols + (["PLPEE01-DVBVD01", "DVHAR01-DVBVD01"] if kind == "detrended" else []):
            for lag in lags:  # positive lag: SOI leads sea level
                r, p, ne, n = corr_sig(a_s[st].loc["1991":], s_s.shift(lag).loc["1991":])
                rows.append(dict(anomaly=kind, station=st, smooth_months=smooth, lag_months=lag,
                                 r=r, p=p, n=n, n_eff=ne))
    C = pd.DataFrame(rows)
    C.round(4).to_csv(os.path.join(OUT, "soi_lag_correlation.csv"), index=False)
    best = C.loc[C.groupby(["anomaly", "station", "smooth_months"])["r"].idxmax()].set_index(
        ["anomaly", "station", "smooth_months"])

    # regression slope (cm per SOI unit), 5-month smoothed at each station's best lag
    a5 = anom_dt.rolling(5, center=True, min_periods=5).mean()
    s5 = soi.rolling(5, center=True, min_periods=5).mean()
    slopes = {}
    for st in cols:
        lag = int(best.loc[("detrended", st, 5), "lag_months"])
        x, y = s5.shift(lag), a5[st]
        m = x.notna() & y.notna() & (x.index >= "1991")
        slopes[st] = dict(lag=lag, cm_per_soi=100 * np.polyfit(x[m], y[m], 1)[0])

    # ---- 2026 status
    recent = {}
    for st in ["HILLARYS", "FFFBH01", "DVBVD01", "MHMAN02", "PLPEE01", "DVCAD01", "DVHAR01"]:
        s = monthly[st].loc["2026"].dropna()
        a = anom_dt[st].loc["2026"].dropna()
        clim_all = monthly[st].loc[CLIM[0]:CLIM[1]]
        pct = {f"{i:%Y-%m}": float((clim_all[clim_all.index.month == i.month] < v).mean() * 100)
               for i, v in s.items()}
        recent[st] = dict(
            last_month=f"{s.index.max():%Y-%m}",
            mean_2026_mAHD=float(s.mean()),
            anom_detrended_cm={f"{i:%Y-%m}": round(float(v) * 100, 1) for i, v in a.items()},
            percentile_vs_1996_2025=pct)
    # lowest Jun-Oct seasonal anomaly ranking
    wint = anom_dt.loc[anom_dt.index.month.isin([6, 7, 8, 9])]
    wint_y = wint.groupby(wint.index.year).mean()
    rank = {st: dict(jun_sep_2026_anom_cm=round(float(wint_y.loc[2026, st]) * 100, 1),
                     rank_lowest=int(wint_y[st].dropna().rank().loc[2026]),
                     of_years=int(wint_y[st].notna().sum()))
            for st in ["HILLARYS", "FFFBH01", "DVBVD01", "MHMAN02", "PLPEE01", "DVCAD01", "DVHAR01",
                       "PLPEE01-DVBVD01", "DVHAR01-DVBVD01"]
            if not np.isnan(wint_y.loc[2026, st])}

    summary = dict(
        datum_steps=[{k: (round(v, 4) if isinstance(v, float) else v) for k, v in s.items()} for s in steps],
        ahd_offsets_m={k: round(float(v), 4) for k, v in off.items()},
        hillarys_minus_fremantle_AHD_cm=round(float(((Mh["HILLARYS"] - 0.763) - (Mh["FFFBH01"] - 0.756)).mean() * 100), 2),
        trend_mm_per_yr_1996_2025={k: round(v, 2) for k, v in trend.items()},
        best_correlation={f"{kind}|{st}|{sm}": dict(lag=int(r.lag_months), r=round(r.r, 3), p=float(f"{r.p:.2g}"),
                                             n_eff=round(r.n_eff, 1))
                          for (kind, st, sm), r in best.iterrows()},
        regression_5mo=slopes, status_2026=recent, jun_sep_2026_rank=rank,
        soi_2026={f"{i:%Y-%m}": float(v) for i, v in soi.loc["2026"].dropna().items()},
    )
    with open(os.path.join(OUT, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    figures(ahd, sub, monthly, anom_dt, soi, C[C.anomaly == "detrended"], clim)
    print(json.dumps(summary, indent=1))


# ----------------------------------------------------------------------------- figures
def label_end(ax, x, y, text, color):
    ax.annotate(text, (x, y), xytext=(4, 0), textcoords="offset points", va="center",
                fontsize=8, color=INK)
    ax.plot([x], [y], "o", ms=4, color=color, mec="white", mew=1)


def figures(ahd, sub, monthly, anom, soi, C, clim):
    # Fig 1: filter demo, Peel Inlet 2026
    fig, ax = plt.subplots(figsize=(9, 3.2))
    w = slice("2026-05-01", "2026-07-31")
    ax.plot(ahd.loc[w, "PLPEE01"].index, ahd.loc[w, "PLPEE01"], color="#b6d5f8", lw=0.8,
            label="Hourly observed")
    ax.plot(sub.loc[w, "PLPEE01"].index, sub.loc[w, "PLPEE01"], color=COL["PLPEE01"], lw=2,
            label="40-h Butterworth low-pass")
    ax.plot(sub.loc[w, "DVBVD01"].index, sub.loc[w, "DVBVD01"], color=COL["DVBVD01"], lw=2,
            label="Cape Bouvard (ocean), low-pass")
    ax.set_ylabel("Water level (m AHD, approx.)")
    ax.set_title("Peel Inlet: tides removed by a 40-hour low-pass filter (May–Jul 2026)",
                 loc="left", color=INK)
    ax.legend(loc="upper left", ncol=3)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
    fig.savefig(os.path.join(FIG, "fig1_filter_example.png"))
    plt.close(fig)

    # Fig 2: monthly subtidal anomalies + SOI, two stacked panels (no dual axis)
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True,
                                 gridspec_kw=dict(height_ratios=[2.2, 1], hspace=0.08))
    a5 = anom.rolling(5, center=True, min_periods=3).mean()
    for st in ["HILLARYS", "DVBVD01", "MHMAN02", "PLPEE01", "DVCAD01", "DVHAR01"]:
        s = a5[st].loc["1991":] * 100
        a1.plot(s.index, s, color=COL[st], lw=1.3, label=NAMES[st])
    a1.axhline(0, color=INK2, lw=0.8)
    a1.set_ylabel("Detrended subtidal anomaly (cm)\n5-month running mean")
    a1.set_title("Peel-Harvey and ocean subtidal water level anomalies vs the Southern Oscillation Index",
                 loc="left", color=INK, pad=26)
    a1.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=6, fontsize=7.5, handlelength=1.5,
              columnspacing=1.2)
    s5 = soi.rolling(5, center=True, min_periods=3).mean().loc["1991":]
    a2.fill_between(s5.index, 0, s5, where=s5 >= 0, color="#257de0", lw=0, label="La Niña-like (SOI > 0)")
    a2.fill_between(s5.index, 0, s5, where=s5 < 0, color="#ec4647", lw=0, label="El Niño-like (SOI < 0)")
    a2.axhline(0, color=INK2, lw=0.8)
    a2.set_ylabel("SOI (BoM)\n5-month mean")
    a2.legend(loc="lower left", ncol=2, fontsize=8)
    a2.xaxis.set_major_locator(mdates.YearLocator(5))
    a2.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    fig.savefig(os.path.join(FIG, "fig2_anomaly_vs_soi.png"))
    plt.close(fig)

    # Fig 3: lag correlation (5-month smoothed)
    fig, ax = plt.subplots(figsize=(7, 3.8))
    c5 = C[C.smooth_months == 5]
    for st in ["HILLARYS", "DVBVD01", "MHMAN02", "PLPEE01", "DVCAD01", "DVHAR01"]:
        d = c5[c5.station == st]
        ax.plot(d.lag_months, d.r, color=COL[st], lw=2, marker="o", ms=4, mec="white", label=NAMES[st])
    ax.axvline(0, color=INK2, lw=0.8)
    ax.set_xlabel("Lag (months; positive = SOI leads water level)")
    ax.set_ylabel("Correlation r")
    ax.set_title("Lagged correlation, SOI vs subtidal level anomaly (5-mo means, 1991–2026)",
                 loc="left", color=INK)
    ax.legend(fontsize=8, loc="lower center", ncol=2)
    fig.savefig(os.path.join(FIG, "fig3_lag_correlation.png"))
    plt.close(fig)

    # Fig 4: 2025-2026 daily subtidal vs climatology band, Peel Inlet & Mandurah & Bouvard
    fig, axes = plt.subplots(3, 1, figsize=(10, 7.5), sharex=True, gridspec_kw=dict(hspace=0.12))
    for ax, st in zip(axes, ["DVBVD01", "MHMAN02", "PLPEE01"]):
        d = sub[st].resample("D").mean()
        base = d.loc[CLIM[0]:CLIM[1]]
        doy = base.groupby(base.index.dayofyear)
        q = pd.DataFrame({"p10": doy.quantile(0.1), "p50": doy.quantile(0.5), "p90": doy.quantile(0.9)})
        q = q.rolling(15, center=True, min_periods=1).mean()
        idx = pd.date_range("2025-01-01", d.index.max(), freq="D")
        qq = q.reindex(idx.dayofyear).set_index(idx)
        ax.fill_between(idx, qq.p10, qq.p90, color="#e6e6e6", lw=0, label="1996–2025 10–90th pct")
        ax.plot(idx, qq.p50, color=INK2, lw=1, ls="--", label="1996–2025 median")
        ax.plot(d.loc[idx[0]:].index, d.loc[idx[0]:], color=COL[st], lw=1.6, label=f"{NAMES[st]} daily subtidal")
        ax.set_ylabel("m AHD (approx.)")
        ax.legend(loc="upper left", ncol=3, fontsize=8)
    axes[0].set_title("2025–2026 subtidal water level against the 1996–2025 seasonal range",
                      loc="left", color=INK)
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    fig.savefig(os.path.join(FIG, "fig4_2026_vs_climatology.png"))
    plt.close(fig)

    # Fig 5: scatter, Peel Inlet anomaly vs SOI (5-mo), 2026 highlighted
    a5 = anom.rolling(5, center=True, min_periods=5).mean()
    s5 = soi.rolling(5, center=True, min_periods=5).mean()
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    m = a5["PLPEE01"].notna() & s5.notna() & (a5.index >= "1995")
    ax.scatter(s5[m], a5["PLPEE01"][m] * 100, s=12, color="#84b9f3", lw=0, label="1995–2025")
    m26 = m & (a5.index >= "2026")
    ax.scatter(s5[m26], a5["PLPEE01"][m26] * 100, s=40, color="#ec4647", edgecolor="white", lw=1, label="2026")
    k = np.polyfit(s5[m], a5["PLPEE01"][m] * 100, 1)
    xx = np.linspace(s5[m].min(), s5[m].max(), 10)
    ax.plot(xx, np.polyval(k, xx), color=INK, lw=1.2)
    ax.annotate(f"{k[0]:.2f} cm per SOI unit", (xx[-1], np.polyval(k, xx[-1])), xytext=(-110, 10),
                textcoords="offset points", fontsize=8, color=INK)
    ax.axhline(0, color=INK2, lw=0.6)
    ax.axvline(0, color=INK2, lw=0.6)
    ax.set_xlabel("SOI, 5-month mean (zero lag)")
    ax.set_ylabel("Peel Inlet subtidal anomaly (cm)")
    ax.set_title("Peel Inlet level vs SOI", loc="left", color=INK)
    ax.legend(loc="upper left", fontsize=8)
    fig.savefig(os.path.join(FIG, "fig5_scatter_peel_soi.png"))
    plt.close(fig)


if __name__ == "__main__":
    main()
