"""Write data/stations.csv: the unique gauge stations in data/gauges.csv, sorted by (src, sid).
Reproduces the July file exactly (its original producer was not saved; logged in T10_RESULTS.md)."""
import pandas as pd

g = pd.read_csv("data/gauges.csv")
g.drop_duplicates(["src", "sid"])[["src", "sid", "lat", "lon"]].sort_values(["src", "sid"]).to_csv("data/stations.csv", index=False)
