"""Loaders for DoT WA tide files, BoM ABSLMP hourly (Hillarys) and UHSLC hourly.

All loaders return a pandas Series of water level in metres, indexed by UTC time.
"""
import glob
import os
import re

import numpy as np
import pandas as pd

ROOT = os.path.join(os.path.dirname(__file__), "..", "data")
DOT = os.path.join(ROOT, "dot")

# DoT LOCATION_ID -> display name (from DoT TIDES_STATIONS_ldot feature service)
STATIONS = {
    "MHMAN02": "Mandurah (Ocean Marina)",
    "MHMAN01": "Mandurah (Fishermans Jetty)",
    "PLPEE01": "Peel Inlet",
    "DVCAD01": "Caddadup",
    "DVHAR01": "Harvey",
    "FFFBH01": "Fremantle FBH",
}

_REC = re.compile(rb"(-?\d+)[ \t]*,[ \t]*(\d{8})\.(\d{4})\d*")


def read_dot_file(path: str) -> pd.Series:
    """Parse one DoT 'H(cm),YYYYMMDD.hhmm/' file. Records may be separated by
    newlines, bare CRs or only '/'. Times are AWST and converted to UTC."""
    raw = open(path, "rb").read()
    body = raw[raw.find(b"YYYYMMDD.hhmm"):]  # skip header
    recs = _REC.findall(body)
    vals = [int(v) for v, _, _ in recs]
    stamps = [d + hm for _, d, hm in recs]
    if not vals:
        return pd.Series(dtype=float)
    t = pd.to_datetime(pd.Series(stamps).str.decode("ascii"), format="%Y%m%d%H%M",
                       errors="coerce")
    ok = t.notna().to_numpy()  # a few files contain invalid stamps (e.g. hour 24)
    s = pd.Series(np.asarray(vals, float)[ok] / 100.0,
                  index=pd.DatetimeIndex(t[ok]) - pd.Timedelta(hours=8))
    s[s <= -99] = np.nan  # -9999 cm flag
    return s


def load_dot(station: str) -> pd.Series:
    """Concatenate a DoT station's historical packet + 2026 recent file(s).

    For stations with paired sensors (files ending a/b, e.g. Fremantle), the
    'a' sensor is used and gaps filled from 'b'.
    """
    files = sorted(glob.glob(os.path.join(DOT, station, f"{station}_0*", "*.txt")))
    files += sorted(glob.glob(os.path.join(DOT, station, "*_2026*.txt")))
    files += sorted(glob.glob(os.path.join(DOT, station, "recent", "*.txt")))
    prim, sec = [], []
    for f in files:
        s = read_dot_file(f)
        base = os.path.splitext(os.path.basename(f))[0]
        (sec if re.search(r"\d{4}b$", base) else prim).append(s)
    p = pd.concat(prim).sort_index()
    p = p[~p.index.duplicated(keep="last")]
    if sec:
        b = pd.concat(sec).sort_index()
        b = b[~b.index.duplicated(keep="last")]
        p = p.combine_first(b)
    p.name = station
    return p


def load_hillarys() -> pd.Series:
    """BoM ABSLMP IDO71012 hourly sea level (m, UTC). -9999 flags removed."""
    out = []
    for f in sorted(glob.glob(os.path.join(ROOT, "hillarys", "IDO71012_*.csv"))):
        df = pd.read_csv(f, skipinitialspace=True)
        df.columns = [c.strip() for c in df.columns]
        t = pd.to_datetime(df.iloc[:, 0].str.strip(), format="%d-%b-%Y %H:%M")
        s = pd.Series(pd.to_numeric(df["Sea Level"], errors="coerce").to_numpy(), index=t)
        out.append(s)
    s = pd.concat(out).sort_index()
    s = s[~s.index.duplicated()]
    s[(s < -50) | (s > 50)] = np.nan
    s.name = "HILLARYS"
    return s


def load_hillarys_pressure() -> pd.Series:
    out = []
    for f in sorted(glob.glob(os.path.join(ROOT, "hillarys", "IDO71012_*.csv"))):
        df = pd.read_csv(f, skipinitialspace=True)
        df.columns = [c.strip() for c in df.columns]
        t = pd.to_datetime(df.iloc[:, 0].str.strip(), format="%d-%b-%Y %H:%M")
        out.append(pd.Series(pd.to_numeric(df["Barometric Pressure"], errors="coerce").to_numpy(), index=t))
    s = pd.concat(out).sort_index()
    s = s[~s.index.duplicated()]
    s[(s < 900) | (s > 1100)] = np.nan
    return s


def load_uhslc(path: str) -> pd.Series:
    """UHSLC hourly CSV: year,month,day,hour,mm (UTC), -32767 missing."""
    df = pd.read_csv(path, header=None, names=["y", "m", "d", "h", "v"])
    t = pd.to_datetime(dict(year=df.y, month=df.m, day=df.d, hour=df.h))
    s = pd.Series(df.v.to_numpy(float) / 1000.0, index=t)
    s[df.v.to_numpy() <= -32767] = np.nan
    return s
