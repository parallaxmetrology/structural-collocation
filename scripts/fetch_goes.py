"""Stream GOES-16 ABI CMIPC band 13 scans (minutes 01/16/31/46), average Tb per hour,
bin to the 0.1-degree MRMS grid (north-first 350 x 700). Raw files never stored.
Output: data/goes/YYYYmmddHH.npz with Tb (K, NaN = <4 pixels) for accumulation hour ending HH.
Plan: PROBE_A2_PLAN.md
"""
import io, os, re, time, urllib.request
from concurrent.futures import ThreadPoolExecutor
import numpy as np, pandas as pd
import h5py
import period as PERIOD

BUCKET = "https://noaa-goes16.s3.amazonaws.com/"
HOURS = PERIOD.HOURS
os.makedirs("data/goes", exist_ok=True)


def get(url):
    for k in range(5):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read()
        except Exception:
            if k == 4:
                raise
            time.sleep(3 * (k + 1))


def grid_index(f):
    """Pixel -> flat 0.1-degree cell index (or -1), via the GOES fixed-grid equations."""
    p = f["goes_imager_projection"].attrs
    H = float(p["perspective_point_height"][0]) + float(p["semi_major_axis"][0])
    a, b = float(p["semi_major_axis"][0]), float(p["semi_minor_axis"][0])
    lon0 = np.deg2rad(float(p["longitude_of_projection_origin"][0]))
    sc = lambda v: v[:].astype(np.float64) * v.attrs["scale_factor"][0] + v.attrs["add_offset"][0]
    x, y = np.meshgrid(sc(f["x"]), sc(f["y"]))
    A = np.sin(x) ** 2 + np.cos(x) ** 2 * (np.cos(y) ** 2 + (a / b) ** 2 * np.sin(y) ** 2)
    B = -2 * H * np.cos(x) * np.cos(y)
    C = H ** 2 - a ** 2
    with np.errstate(invalid="ignore"):
        rs = (-B - np.sqrt(B ** 2 - 4 * A * C)) / (2 * A)
        sx = rs * np.cos(x) * np.cos(y); sy = -rs * np.sin(x); sz = rs * np.cos(x) * np.sin(y)
        lat = np.degrees(np.arctan((a / b) ** 2 * sz / np.sqrt((H - sx) ** 2 + sy ** 2)))
        lon = np.degrees(lon0 - np.arctan(sy / (H - sx)))
    r = np.floor((55.0 - lat) / 0.1); c = np.floor((lon + 130.0) / 0.1)
    ok = np.isfinite(r) & np.isfinite(c) & (r >= 0) & (r < 350) & (c >= 0) & (c < 700)
    idx = np.full(lat.shape, -1, dtype=np.int64)
    idx[ok] = (r[ok] * 700 + c[ok]).astype(np.int64)
    return idx.ravel()


IDX = None


def tb(b):
    with h5py.File(io.BytesIO(b), "r") as f:
        global IDX
        if IDX is None:
            IDX = grid_index(f)
        v = f["CMI"]
        raw = v[:]
        out = raw.astype(np.float32) * v.attrs["scale_factor"][0] + v.attrs["add_offset"][0]
        out[raw == v.attrs["_FillValue"][0]] = np.nan
    return out.ravel()


def one(t):
    out = f"data/goes/{t:%Y%m%d%H}.npz"
    if os.path.exists(out):
        return t, "skip"
    t0 = t - pd.Timedelta(hours=1)
    lst = get(f"{BUCKET}?list-type=2&prefix=ABI-L2-CMIPC/{t0:%Y}/{t0.dayofyear:03d}/{t0:%H}/OR_ABI-L2-CMIPC-M6C13").decode()
    keys = re.findall(r"<Key>([^<]*)</Key>", lst)
    pick = [k for k in keys if re.search(r"_s\d{9}(01|16|31|46)", k)]
    scans = [tb(get(BUCKET + k)) for k in pick]
    if len(scans) < 3:
        return t, f"only {len(scans)} scans"
    v = np.nanmean(np.stack(scans), axis=0)
    ok = (IDX >= 0) & np.isfinite(v)
    n = np.bincount(IDX[ok], minlength=350 * 700)
    s = np.bincount(IDX[ok], v[ok], minlength=350 * 700)
    with np.errstate(invalid="ignore", divide="ignore"):
        cell = np.where(n >= 4, s / n, np.nan).reshape(350, 700).astype(np.float32)
    np.savez_compressed(out, Tb=cell, nscan=len(scans))
    return t, "ok"


if __name__ == "__main__":
    # initialise the pixel->cell map once before threading
    one(HOURS[0])
    n = 0
    with ThreadPoolExecutor(6) as ex:
        for t, s in ex.map(one, HOURS):
            n += 1
            if s != "ok" and s != "skip" or n % 48 == 0:
                print(n, t, s, flush=True)
    print("done", n)
