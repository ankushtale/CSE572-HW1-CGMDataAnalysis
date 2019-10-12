from datetime import datetime, timedelta
import pandas as pd
import pandas


cgm_timestamps = pd.read_csv("dm_proj1_data/CGMDatenumLunchPat1.csv")
cgm_values = pd.read_csv("dm_proj1_data/CGMSeriesLunchPat1.csv")


def mdate_to_pdate(mdate):
    return (datetime.fromordinal(int(mdate)) + timedelta(days=mdate % 1) - timedelta(days=366)).strftime("%Y-%m-%d %H:%M:%S")


for col in cgm_timestamps.columns:
    cgm_timestamps[col] = cgm_timestamps[col].apply(lambda x: mdate_to_pdate(x) if pd.notnull(x) else x)

cgm_timestamps, cgm_values = cgm_timestamps.T, cgm_values.T


cgm_timestamps = cgm_timestamps.set_index([pandas.Index(['cgm_{}'.format(i) for i in range(1, len(cgm_timestamps) + 1)])])
cgm_values = cgm_values.set_index([pandas.Index(['cgm_{}'.format(i) for i in range(1, len(cgm_timestamps) + 1)])])

cgm_merged = pandas.merge(cgm_timestamps, cgm_values, how="inner", left_index=True, right_index=True)

cols = cgm_merged.columns.tolist()
new_cols = []
for i in range(len(cgm_merged)):
    new_cols.append(cols[0 + i])
    new_cols.append(cols[len(cgm_merged.columns) // 2 + i])

cgm_merged = cgm_merged[new_cols]
cgm_merged = cgm_merged.fillna(method='ffill', limit=2)

# sort timeseries
for i in cgm_merged.columns:
    cgm_merged[i] = cgm_merged[i].values[::-1]
