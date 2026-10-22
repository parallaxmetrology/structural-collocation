"""Stream MRMS RadarOnly_QPE_01H and RadarAccumulationQualityIndex_01H (top of hour),
reduce in memory, never store raw GRIB.

Per hour -> data/mrms/YYYYmmddHH.npz with
  R, Q : 0.1-degree block means (350 x 700, north-first), NaN if <80% valid pixels
  Rg, Qg : values at each gauge's own 0.01-degree pixel (order of data/stations.csv)
"""
import gzip, os, sys, urllib.request, time
from concurrent.futures import ThreadPoolExecutor
import numpy as np, pandas as pd
import eccodes
import period as PERIOD

BASE = "https://noaa-mrms-pds.s3.amazonaws.com/CONUS/"
PROD = {"R": "RadarOnly_QPE_01H_00.00", "Q": "RadarAccumulationQualityIndex_01H_00.00"}
HOURS = PERIOD.HOURS
os.makedirs("data/mrms", exist_ok=True)

st = pd.read_csv("data/stations.csv")
pr = np.clip(np.round((54.995 - st.lat.values) / 0.01).astype(int), 0, 3499)
pc = np.clip(np.round((st.lon.values + 129.995) / 0.01).astype(int), 0, 6999)


def grab(prod, t):
    url = f"{BASE}{prod}/{t:%Y%m%d}/MRMS_{prod}_{t:%Y%m%d-%H%M%S}.grib2.gz"
    for k in range(5):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                b = gzip.decompress(r.read())
            break
        except Exception:
            if k == 4:
                raise
            time.sleep(3 * (k + 1))
    h = eccodes.codes_new_from_message(b)
    v = eccodes.codes_get_values(h).reshape(3500, 7000).astype(np.float32)
    eccodes.codes_release(h)
    v[v < 0] = np.nan  # -3 (radar) / -1 (RQI) = no coverage
    return v


def block(v):
    b = v.reshape(350, 10, 700, 10)
    valid = np.isfinite(b).sum(axis=(1, 3))
    m = np.nanmean(b, axis=(1, 3))
    m[valid < 80] = np.nan
    return m.astype(np.float32)


def one(t):
    out = f"data/mrms/{t:%Y%m%d%H}.npz"
    if os.path.exists(out):
        return t, "skip"
    R = grab(PROD["R"], t); Q = grab(PROD["Q"], t)
    with np.errstate(all="ignore"):
        np.savez_compressed(out, R=block(R), Q=block(Q), Rg=R[pr, pc], Qg=Q[pr, pc])
    return t, "ok"


if __name__ == "__main__":
    n = 0
    with ThreadPoolExecutor(6) as ex:
        for t, s in ex.map(one, HOURS):
            n += 1
            if n % 24 == 0:
                print(n, t, s, flush=True)
    print("done", n)
