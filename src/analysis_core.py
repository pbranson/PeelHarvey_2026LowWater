"""Core routines: 40-h Butterworth low-pass (subtidal) filter and SOI loading."""
import numpy as np
import pandas as pd
from scipy import signal


def to_hourly(series: pd.Series, max_gap_h: int = 6) -> pd.Series:
    """Resample to an hourly grid (hour-centred mean) and interpolate short gaps only."""
    s = series.sort_index()
    s = s[~s.index.duplicated()]
    h = s.resample("1h").mean()
    return h.interpolate(limit=max_gap_h, limit_area="inside")


def butter_lowpass(series: pd.Series, cutoff_h: float = 40.0, order: int = 4,
                   min_segment_h: int = 24 * 10) -> pd.Series:
    """Zero-phase (forward-backward) Butterworth low-pass on an hourly series.

    Applied separately to each contiguous (gap-free) segment; segments shorter
    than ``min_segment_h`` are dropped, and a cutoff-length edge zone at each
    segment end is masked to suppress filter start-up transients.
    """
    dt_h = (series.index[1] - series.index[0]) / pd.Timedelta("1h")
    fs = 1.0 / dt_h                     # samples per hour
    wn = (1.0 / cutoff_h) / (fs / 2.0)  # normalised cutoff
    sos = signal.butter(order, wn, btype="low", output="sos")
    out = pd.Series(np.nan, index=series.index)
    valid = series.notna().to_numpy()
    edges = np.flatnonzero(np.diff(np.r_[0, valid.astype(int), 0]))
    trim = int(round(cutoff_h / dt_h))
    for a, b in zip(edges[::2], edges[1::2]):
        if b - a < min_segment_h:
            continue
        y = signal.sosfiltfilt(sos, series.iloc[a:b].to_numpy())
        y[:trim] = np.nan
        y[-trim:] = np.nan
        out.iloc[a:b] = y
    return out


def load_soi_bom(path: str) -> pd.Series:
    """BoM monthly SOI (IDCKGSM000, 'YYYYMM,value' lines) -> Series on month start."""
    df = pd.read_csv(path, header=None, names=["ym", "soi"], dtype={"ym": str})
    df["ym"] = df["ym"].str.strip()
    df["date"] = pd.to_datetime(df["ym"], format="%Y%m")
    return df.set_index("date")["soi"].astype(float)
