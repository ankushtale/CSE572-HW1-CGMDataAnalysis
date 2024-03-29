from datetime import datetime, timedelta
import pandas as pd
import pandas
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from scipy.fftpack import fft

numFiles = 5
cgm_values = [pd.read_csv('dm_proj1_data/CGMSeriesLunchPat{}.csv'.format(i)) for i in range(1, numFiles + 1)]
cgm_timestamps = [pd.read_csv('dm_proj1_data/CGMDatenumLunchPat{}.csv'.format(i)) for i in range(1, numFiles + 1)]


def mdate_to_pdate(mdate):
    return (datetime.fromordinal(int(mdate)) + timedelta(days=mdate % 1) - timedelta(days=366)).strftime("%Y-%m-%d %H:%M:%S")


def epochToDate(cgmT):
    for col in cgmT.columns:
        cgmT[col] = cgmT[col].apply(lambda x: mdate_to_pdate(x) if pd.notnull(x) else x)


for i in range(numFiles):
    epochToDate(cgm_timestamps[i])

for i in range(numFiles):
    cgm_timestamps[i], cgm_values[i] = cgm_timestamps[i].T, cgm_values[i].T

for i in range(numFiles):
    cgm_timestamps[i] = cgm_timestamps[i].set_index([pandas.Index(['cgm_{}'.format(i) for i in range(1, len(cgm_timestamps[i]) + 1)])])
    cgm_values[i] = cgm_values[i].set_index([pandas.Index(['cgm_{}'.format(i) for i in range(1, len(cgm_timestamps[i]) + 1)])])


cgm_merged = []
for i in range(numFiles):
    cgm_merged.append(pandas.merge(cgm_timestamps[i], cgm_values[i], how="inner", left_index=True, right_index=True))


def colNames(cgmM):
    cols = cgmM.columns.tolist()
    new_cols = []
    for i in range(len(cols) // 2):
        new_cols.extend([cols[i], cols[len(cols) // 2 + i]])
    return new_cols


cols = []
for i in range(numFiles):
    cols.append(colNames(cgm_merged[i]))


for i in range(numFiles):
    cgm_merged[i] = cgm_merged[i][cols[i]]
    cgm_merged[i] = cgm_merged[i].fillna(method='ffill', limit=2)


# sort timeseries
for j in range(numFiles):
    for i in cgm_merged[j].columns:
        cgm_merged[j][i] = cgm_merged[j][i].values[::-1]


def addOutofRange(cgmMO):
    num = 2
    for i in range(len(cgmMO.columns) // 2):
        cgmMO.insert(num, '{}_out_of_range'.format(i), [70 - j if j < 70 else j - 180 if j > 180 else 0 for j in cgmMO['{}_y'.format(i)]], True)
        num += 3


def delOutOfRange(cgmMO):
    out = []
    for i in range(len(cgmMO.columns) // 3):
        out.append(round(sum(cgmMO['{}_out_of_range'.format(i)]) / len(cgmMO['{}_out_of_range'.format(i)]), 2))
        del cgmMO['{}_out_of_range'.format(i)]
    return out


for i in range(numFiles):
    addOutofRange(cgm_merged[i])


def countOutOfRange(cgmMO):
    posCount, negCount = 0, 0
    outOfRange = []
    for i in range(len(cgmMO.columns) // 3):
        for j in cgmMO['{}_out_of_range'.format(i)]:
            if j > 0:
                posCount += 1
            if j < 0:
                negCount += 1
    #     outOfRange.append([posCount,negCount])
        outOfRange.append(posCount)
        posCount, negCount = 0, 0
    return outOfRange


outOfRangeList = []
for i in range(numFiles):
    outOfRangeList.append([countOutOfRange(cgm_merged[i]), delOutOfRange(cgm_merged[i])])


timeRange = [
    [
        '{1} - {0}'.format(cgm_merged[idx]['{}_x'.format(i)][0], cgm_merged[idx]['{}_x'.format(i)][-1]) for i in range(len(cgm_merged[idx].columns) // 2)
    ] for idx in range(numFiles)
]

feature_matrix = [pd.DataFrame({'timeRange': timeRange[i], '#OutOfRange': outOfRangeList[i][0], 'meanOutOfRange':outOfRangeList[i][1]}) for i in range(numFiles)]
