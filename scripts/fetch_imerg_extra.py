"""Second pass over the same IMERG Late files: keep the separate infrared and microwave
estimates and IMERG's own error field (cropped). Output: data/imerg_x/YYYYmmddHHMM.npz (350 x 700,
north-first, mm/h, NaN = fill). Token read from secrets/earthdata_token.txt, never printed.
"""
import io, os, re, time, urllib.request
from concurrent.futures import ThreadPoolExecutor
import numpy as np, pandas as pd, h5py
import period as PERIOD

TOK = "".join(l.strip() for l in open("secrets/earthdata_token.txt") if not l.startswith("#"))
BASE = "https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L3/GPM_3IMERGHHL.07/2024/"
os.makedirs("data/imerg_x", exist_ok=True)


def get(url, auth=True):
    for k in range(6):
        try:
            req = urllib.request.Request(url, headers={"Authorization": f"Bearer {TOK}"} if auth else {})
            with urllib.request.urlopen(req, timeout=180) as r:
                return r.read()
        except Exception as e:
            if k == 5:
                raise RuntimeError(f"failed {url.split('/')[-1]}: {type(e).__name__}") from None
            time.sleep(10 * (k + 1))


def one(name):
    s = re.search(r"\.(\d{8})-S(\d{4})", name)
    out = f"data/imerg_x/{s.group(1)}{s.group(2)}.npz"
    if os.path.exists(out):
        return name, "skip"
    doy = pd.Timestamp(s.group(1)).dayofyear
    b = get(f"{BASE}{doy:03d}/{name}")
    with h5py.File(io.BytesIO(b), "r") as f:
        lat = f["Grid/lat"][1100:1450]; lon = f["Grid/lon"][500:1200]
        assert abs(lat[0] - 20.05) < 1e-3 and abs(lon[0] + 129.95) < 1e-3, (lat[0], lon[0])
        def crop(path):
            a = f[path][0, 500:1200, 1100:1450].astype(np.float32)
            a[a < 0] = np.nan
            return a.T[::-1]
        np.savez_compressed(out, IR=crop("Grid/Intermediate/IRprecipitation"),
                            MW=crop("Grid/Intermediate/MWprecipitation"),
                            ERR=crop("Grid/randomError"))
    return name, "ok"


if __name__ == "__main__":
    names = []
    for day in PERIOD.DAYS:
        html = get(f"{BASE}{day.dayofyear:03d}/", auth=False).decode()
        names += sorted(set(re.findall(r'(3B-HHR-L\.MS\.MRG\.3IMERG\.\d{8}-S\d{6}-E\d{6}\.\d{4}\.V07B\.HDF5)"', html)))
    print("files", len(names), flush=True)
    n = 0
    with ThreadPoolExecutor(2) as ex:
        for name, st in ex.map(one, names):
            n += 1
            if n % 48 == 0:
                print(n, name[22:40], st, flush=True)
    print("done", n)
