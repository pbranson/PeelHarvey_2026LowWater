# Peel-Harvey Estuary subtidal water levels and the SOI

Analysis of low-frequency (subtidal) water level in the Peel-Harvey Estuary, WA,
and its relationship with the Southern Oscillation Index (SOI), the ENSO index
that modulates the strength of the Leeuwin Current and hence coastal sea level
along south-west WA.

## Data

| Gauge | Source | ID | Record used | Datum in file |
|---|---|---|---|---|
| Peel Inlet | DoT WA | PLPEE01 | 1985 – Oct 2026 | CD (2017+, 1.201 m below BM A895); earlier "not available" |
| Caddadup | DoT WA | DVCAD01 | 1994 – Sep 2026 | "Under review" |
| Harvey | DoT WA | DVHAR01 | 1989 – Sep 2026 | "Not available" |
| Mandurah Ocean Marina | DoT WA | MHMAN02 | 2007 – Sep 2026 | CD, 2.682 m below BM MAN TIDAL 02 |
| Mandurah Fishermans Jetty | DoT WA | MHMAN01 | 1990 – 2007 (spliced before MHMAN02) | CD, 2.065 m below TBM 1971 |
| Cape Bouvard (open coast) | DoT WA | DVBVD01 | 1995 – Sep 2026 | CD, 16.251 m below BM Florida 2 |
| Fremantle FBH | DoT WA | FFFBH01 | 1986 – Oct 2026 | LWM Fremantle 1949 (= AHD − 0.756 m) |
| Hillarys | BoM ABSLMP | IDO71012 | 1991 – Sep 2026 (hourly) | Gauge zero; AHD = 0.763 m |
| SOI | BoM | IDCKGSM000 | monthly to Oct 2026 | – |

The DoT files come from the S3 links behind DoT's *Download tide and wave data*
map (the `TIDES_STATIONS_ldot` ArcGIS feature service). The 2026 "recent data"
files cover 1 Jan to 30 Sep (Peel Inlet to 31 Oct) at 5-minute resolution.
UHSLC Fremantle (h175) and NOAA SOI are also downloaded as cross-checks.
`src/download_data.sh` re-downloads everything.

AODN/IMOS: there is no DoT water-level layer on the AODN geoserver. The relevant
hourly national product is BoM's ABSLMP, which is used here (Hillarys).

## Method

1. **Parsing**: DoT `H(cm),YYYYMMDD.hhmm/` files (mixed CR/LF/`/` record
   separators, −9999 flags, AWST → UTC). For Fremantle's paired sensors (a/b),
   sensor *a* is used and its gaps are filled from *b*.
2. **Hourly grid**: hour means; gaps ≤ 6 h are linearly interpolated.
3. **Datums**, worked out on monthly means (`src/datums.py`):
   - **Reference**: Hillarys (ABSLMP, AHD). The independent check against Fremantle
     LWMF−0.756 m agrees to **1.7 cm** on average (drift ≈ +1 mm/yr, consistent
     with differential land subsidence). This confirms the AHD reductions.
   - **Steps** (station − Hillarys in windows either side of a break). Each record
     is aligned to its latest segment.

     | Station | Break | Step (m) | Note |
     |---|---|---|---|
     | Peel Inlet | Jan 1995 | −0.187 | Estimated on Nov–Apr months to avoid the Dawesville Cut effect |
     | Harvey | Jan 1997 | +0.037 | |
     | Cape Bouvard | 2013–2015 outage | −0.051 | 1993–94 dropped (sensor settling) |
     | Fremantle DoT | Jan 2013 | −0.030 | |

   - **AHD reduction** for stations with no published AHD relation:
     - Open-coast gauges are matched to Hillarys over 2015–2025.
     - Estuary gauges are matched to Mandurah in Jan–Mar (low runoff).

     The estimated offsets (m to add to CD) are:

     | Mandurah | Peel | Caddadup | Harvey | Cape Bouvard |
     |---|---|---|---|---|
     | −0.499 | −0.514 | −0.562 | −0.603 | −0.608 |

     Absolute AHD values for the Peel stations are therefore approximate. Anomalies
     and correlations do not depend on these constants. They can be replaced with
     surveyed benchmark AHD heights (A895, MAN TIDAL 02, Florida 2) when available.
4. **Subtidal filter**: 4th-order Butterworth low-pass, 40-h cutoff, applied
   forward–backward (zero phase; effective 8th order) on each gap-free segment;
   one cutoff length trimmed at segment ends. Amplitude response: M2 0.08 %,
   K1 1.6 %, O1 2.9 %, 40 h 50 %, 3 days 99 %.
5. **Anomalies**: monthly means of the subtidal series (≥ 20 days of data) minus
   the 1996–2025 monthly climatology. The baseline is post-Dawesville Cut, which
   opened in April 1994. A **linear trend** fitted over 1996–2025 is also removed,
   giving the "detrended" anomalies used for the main results.
   Estuary − Cape Bouvard anomalies isolate estuary-specific signals.
6. **SOI correlation**: lags −3 to +12 months (positive = SOI leads), raw and
   5-month running means. p-values use the Bretherton et al. (1999) effective
   sample size. An inverse-barometer-corrected Hillarys series uses Hillarys
   pressure.

## Results (to Sep/Oct 2026)

**Correlation with SOI.** All gauges respond coherently to ENSO. The SOI leads by
about **2 months** everywhere (detrended anomalies, 1991–2026):

| Gauge | r (monthly, lag 2) | r (5-mo mean, best lag) | Slope (cm per SOI unit) |
|---|---|---|---|
| Hillarys | 0.59 | 0.76 | 0.54 |
| Fremantle | 0.59 | 0.77 | 0.55 |
| Cape Bouvard | 0.54 | 0.71 | 0.53 |
| Mandurah | 0.57 | 0.72 | 0.57 |
| Peel Inlet | 0.54 | 0.73 (lag 1) | 0.49 |
| Caddadup | 0.60 | 0.79 | 0.60 |
| Harvey | 0.55 | 0.74 (lag 1) | 0.53 |

All are significant at p < 0.005 for 5-month means and p < 1e−10 for monthly
values; exact values are in `outputs/summary.json`. Estuary − ocean anomalies are
*not* significantly correlated with SOI (r ≤ 0.21, p ≥ 0.07). The ENSO signal
reaches the estuary through the ocean (Leeuwin Current / coastal sea level),
transmitted via the Mandurah Channel and Dawesville Cut.

**2026.** The SOI swung from about +10 (Jan–Mar) to −11 (Apr), −25 (Jun), −29 (Jul)
and −21 (Sep). Detrended subtidal anomalies fell from near zero in May to:

| Gauge | Jul 2026 | Sep 2026 | Jun–Sep mean | Rank (lowest) |
|---|---|---|---|---|
| Peel Inlet | −9.7 cm | −10.1 cm | −7.7 cm | 3rd of 39 years |
| Harvey | −10.7 cm | −11.2 cm | −8.8 cm | **lowest of 35** |
| Caddadup | −9.1 cm | −9.7 cm | −7.6 cm | 4th of 32 |
| Mandurah | −8.7 cm | −10.0 cm | −7.2 cm | 3rd of 36 |
| Fremantle | −9.6 cm | −10.1 cm | −7.9 cm | 2nd of 40 |

The magnitude matches the regression: about 0.55 cm per SOI unit × −20 gives
about −11 cm.

Peel Inlet was also −8 cm in October. Estuary − ocean anomalies added only −1 to
−2 cm in winter 2026, consistent with lower river inflow in an El Niño winter.

Without detrending, the 2026 anomalies are about 4–6 cm smaller. Sea level has
risen about 3–6 mm/yr since 1996 (trend estimates: Fremantle/Hillarys ≈ 5.6–5.7,
Peel stations 2.3–3.6 mm/yr; the difference likely reflects land motion and
residual datum drift). So 2026 levels are low against the trend but less
extreme against the raw 30-year mean.

## Caveats

- The datum steps are statistical estimates, not survey ties. AHD for the Peel
  stations is approximate (± a few cm).
- DoT is moving Peel-area gauges from CD to LAT. The 2026 files are still on CD.
  Future files must be checked for a step.
- Linear detrending over 1996–2025 partly aliases decadal ENSO variability (the
  La Niña-dominated 2010s).
- 5-month-mean significance uses few effective degrees of freedom (n_eff ≈ 12–22).

## Outputs

- `outputs/subtidal_daily_mAHD.csv`: daily subtidal level, all gauges
- `outputs/subtidal_monthly_mAHD.csv`, `subtidal_monthly_anomaly_m.csv`,
  `subtidal_monthly_anomaly_detrended_m.csv`
- `outputs/soi_lag_correlation.csv`: all lags, raw/detrended, 1- and 5-month
- `outputs/summary.json`: datum steps, offsets, trends, best correlations, 2026 status
- `figures/fig1–fig5`

## Reproduce

```
pip install -r requirements.txt
bash src/download_data.sh
python src/run_analysis.py
```
