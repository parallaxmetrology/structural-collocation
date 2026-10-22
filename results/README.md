# Saved results

The analysis outputs as recorded, for July 2024 (`july/`) and the August 2024 replication (`august/`):
result JSONs, the station table, the point-support table (`t7_points.csv`) and the Probe A2 window table.

The scripts write to `data/` (July) and `data_aug/` (August). To rebuild `PAPER.html` from these saved
results without re-fetching inputs:

```
mkdir -p data data_aug
cp results/july/* data/ && cp results/august/* data_aug/
python build_paper.py
```

`build_labbook.py` also hashes the reduced inputs for its digest table, so it needs the fetched data.
