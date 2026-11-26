#!/usr/bin/env bash
# Download all raw data used by the analysis into data/. Re-run to refresh.
set -euo pipefail
cd "$(dirname "$0")/.."
UA="Mozilla/5.0"
S3=https://s3-ap-southeast-2.amazonaws.com/transport.wa/Tide_Packet
S3R=https://s3.ap-southeast-2.amazonaws.com/transport.wa/Tide_Packet/RecentData

# DoT WA station list (ArcGIS feature service behind the DoT download map)
mkdir -p data/dot
curl -sS -o data/dot/tide_stations.json \
  "https://services6.arcgis.com/67Ks15nDmWoIbK8b/arcgis/rest/services/TIDES_STATIONS_ldot/FeatureServer/0/query?where=1%3D1&outFields=*&returnGeometry=true&outSR=4326&f=json"

# DoT WA tide packets (history to end of previous year) + current-year "recent data"
for id in PLPEE01 DVHAR01 DVCAD01 MHMAN02 MHMAN01 DVBVD01 FFFBH01; do
  mkdir -p data/dot/$id
  curl -sS -f -o data/dot/$id/${id}_Tides.zip $S3/${id}_Tides.zip
  (cd data/dot/$id && unzip -o -q ${id}_Tides.zip || true)
done
# Recent-data file names come from the station list (RECENT_DATA_LINK)
python3 - <<'PY'
import json, os, subprocess, zipfile
d = json.load(open("data/dot/tide_stations.json"))
want = {"PLPEE01","DVHAR01","DVCAD01","MHMAN02","DVBVD01","FFFBH01"}
for f in d["features"]:
    a = f["attributes"]; u = a.get("RECENT_DATA_LINK"); sid = a["LOCATION_ID"]
    if sid in want and u:
        dst = os.path.join("data/dot", sid, os.path.basename(u))
        subprocess.run(["curl", "-sS", "-f", "-o", dst, u], check=True)
        if dst.endswith(".zip"):
            zipfile.ZipFile(dst).extractall(os.path.join("data/dot", sid, "recent"))
PY

# BoM ABSLMP hourly, Hillarys (IDO71012)
mkdir -p data/hillarys
for y in $(seq 1991 "$(date +%Y)"); do
  curl -s -f -A "$UA" -o data/hillarys/IDO71012_$y.csv https://www.bom.gov.au/ntc/IDO71012/IDO71012_$y.csv \
    || rm -f data/hillarys/IDO71012_$y.csv
done

# UHSLC Fremantle (station 175) hourly fast-delivery + PSMSL monthly (cross-checks)
mkdir -p data/fremantle
curl -sS -f -o data/fremantle/uhslc_h175.csv https://uhslc.soest.hawaii.edu/data/csv/fast/hourly/h175.csv
curl -sS -f -o data/fremantle/psmsl_111.rlrdata https://psmsl.org/data/obtaining/rlr.monthly.data/111.rlrdata

# SOI: BoM monthly (primary) and NOAA CPC (cross-check)
mkdir -p data/soi
curl -sS -f -A "$UA" -o data/soi/soi_bom_monthly.txt https://www.bom.gov.au/clim_data/IDCKGSM000/soi_monthly.txt
curl -sS -f -o data/soi/soi_noaa.txt https://www.cpc.ncep.noaa.gov/data/indices/soi
rm -f data/processed/hourly_raw.pkl   # force rebuild of cache
echo "done"
