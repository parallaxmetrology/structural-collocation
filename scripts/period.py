"""Period switch (T10_PLAN.md). PRECIP_PERIOD=jul (default) or aug selects the hours and days.
The data root is the working directory's data/: July runs from the project root, August from
runs/aug/ (data -> ../../data_aug, secrets -> ../../secrets)."""
import os
import pandas as pd

SPEC = {"jul": ("2024-07-01 01:00", "2024-07-22 00:00"), "aug": ("2024-08-01 01:00", "2024-08-22 00:00")}
NAME = os.environ.get("PRECIP_PERIOD", "jul")
T0, T1 = (pd.Timestamp(s) for s in SPEC[NAME])
HOURS = pd.date_range(T0, T1, freq="h")
DAYS = pd.date_range(T0.floor("D"), (T1 - pd.Timedelta(hours=1)).floor("D"))
