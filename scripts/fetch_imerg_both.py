"""Fetch efficiency only (T10, logged): download each IMERG Late half-hourly file once and write both
data/imerg (fetch_imerg.one) and data/imerg_x (fetch_imerg_extra.one) with their unchanged crop code."""
import threading
from concurrent.futures import ThreadPoolExecutor
import fetch_imerg as A
import fetch_imerg_extra as B
import period as PERIOD
import re

_orig, _cache, _lock = A.get, {}, threading.Lock()


def cget(url, auth=True):
    if not url.endswith(".HDF5"):
        return _orig(url, auth)
    with _lock:
        if url in _cache:
            return _cache[url]
    b = _orig(url, auth)
    with _lock:
        _cache[url] = b
    return b


A.get = B.get = cget


def both(name):
    try:
        return name, A.one(name)[1] + "/" + B.one(name)[1]
    finally:
        with _lock:
            for k in [k for k in _cache if k.endswith(name)]:
                del _cache[k]


if __name__ == "__main__":
    names = []
    for day in PERIOD.DAYS:
        html = _orig(f"{A.BASE}{day.dayofyear:03d}/", auth=False).decode()
        names += sorted(set(re.findall(r'(3B-HHR-L\.MS\.MRG\.3IMERG\.\d{8}-S\d{6}-E\d{6}\.\d{4}\.V07B\.HDF5)"', html)))
    print("files", len(names), flush=True)
    n = 0
    with ThreadPoolExecutor(4) as ex:
        for name, st in ex.map(both, names):
            n += 1
            if n % 48 == 0:
                print(n, name[22:40], st, flush=True)
    print("done", n)
