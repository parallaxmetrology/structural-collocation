"""Fetch label gauges for the feasibility window: USCRN hourly + ASOS hourly (IEM).

Output: data/gauges.csv  columns: src, sid, lat, lon, t_end (UTC hour ending, ISO), p_mm
Plan: FEASIBILITY_PLAN.md (USCRN P_CALC; ASOS p01i at minute 50-59 -> next hour ending).
"""
import io, re, json, time, urllib.request
import pandas as pd

OUT = "data/gauges.csv"
import period as PERIOD
T0, T1 = PERIOD.T0, PERIOD.T1
LAT, LON = (20.0, 55.0), (-130.0, -60.0)
CONUS = ("AL AZ AR CA CO CT DE FL GA ID IL IN IA KS KY LA ME MD MA MI MN MS MO MT NE NV NH NJ "
         "NM NY NC ND OH OK OR PA RI SC SD TN TX UT VT VA WA WV WI WY DC").split()


def get(url, tries=4):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read()
        except Exception as e:
            if k == tries - 1:
                raise
            time.sleep(5 * (k + 1))


def uscrn():
    base = "https://www.ncei.noaa.gov/pub/data/uscrn/products/hourly02/2024/"
    names = sorted(set(re.findall(r'CRNH0203-2024-([A-Z]{2})_[^"]*\.txt', get(base).decode())))
    files = sorted(set(re.findall(r'(CRNH0203-2024-[A-Z]{2}_[^"]*\.txt)', get(base).decode())))
    rows = []
    for f in files:
        if f[14:16] not in CONUS:
            continue
        txt = get(base + f).decode()
        for line in txt.splitlines():
            p = line.split()
            if len(p) < 13:
                continue
            t = pd.Timestamp(p[1] + p[2].zfill(4))
            if not (T0 <= t <= T1):
                continue
            v = float(p[12])
            if v <= -9999:
                continue
            rows.append(("USCRN", p[0], float(p[7]), float(p[6]), t, v))
    return pd.DataFrame(rows, columns=["src", "sid", "lat", "lon", "t_end", "p_mm"])


def asos():
    frames = []
    for st in CONUS:
        gj = json.loads(get(f"https://mesonet.agron.iastate.edu/geojson/network/{st}_ASOS.geojson"))
        ids = []
        for ft in gj["features"]:
            a = ft["properties"].get("attributes") or {}
            rm = a.get("METAR_RESET_MINUTE")
            if a.get("HAS_PHOUR") == "1" and rm is not None and 50 <= int(rm) <= 59:
                ids.append(ft["id"])
        for i in range(0, len(ids), 40):
            q = "&".join(f"station={s}" for s in ids[i:i + 40])
            url = ("https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?" + q +
                   f"&data=p01i&year1={T0.year}&month1={T0.month}&day1={T0.day}&year2={T1.year}&month2={T1.month}&day2={T1.day}"
                   "&tz=Etc/UTC&format=onlycomma&latlon=yes&missing=M&trace=0.0001&report_type=3")
            df = pd.read_csv(io.BytesIO(get(url)))
            frames.append(df)
            time.sleep(1)
    d = pd.concat(frames)
    d = d[d.p01i != "M"].copy()
    d["valid"] = pd.to_datetime(d["valid"])
    d = d[(d.valid.dt.minute >= 50) & (d.valid.dt.minute <= 59)]
    d["t_end"] = d.valid.dt.floor("h") + pd.Timedelta(hours=1)
    d["p_mm"] = d.p01i.astype(float) * 25.4
    d.loc[d.p_mm < 0.01, "p_mm"] = 0.0  # trace -> 0
    d = d.sort_values("valid").drop_duplicates(["station", "t_end"], keep="last")
    d = d[(d.t_end >= T0) & (d.t_end <= T1)]
    return pd.DataFrame({"src": "ASOS", "sid": d.station, "lat": d.lat, "lon": d.lon,
                         "t_end": d.t_end, "p_mm": d.p_mm})


if __name__ == "__main__":
    u = uscrn(); print("USCRN rows", len(u), "stations", u.sid.nunique())
    a = asos(); print("ASOS rows", len(a), "stations", a.sid.nunique())
    g = pd.concat([u, a])
    g = g[g.lat.between(*LAT) & g.lon.between(*LON)]
    g.to_csv(OUT, index=False)
    print("wrote", OUT, len(g), "rows;", (g.p_mm > 0).mean().round(3), "share wet")
