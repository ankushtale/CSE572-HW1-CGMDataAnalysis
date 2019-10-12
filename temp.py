from datetime import datetime, timedelta
import pandas as pd
import pandas
import matplotlib.pyplot as plt


cgm_timestamps = pd.read_csv("dm_proj1_data/CGMDatenumLunchPat1.csv")
cgm_values = pd.read_csv("dm_proj1_data/CGMSeriesLunchPat1.csv")


def mdate_to_pdate(mdate):
    return (datetime.fromordinal(int(mdate)) + timedelta(days=mdate % 1) - timedelta(days=366)).strftime("%Y-%m-%d %H:%M:%S")


print('first', mdate_to_pdate(737225.5842))
print('last', mdate_to_pdate(737225.48))

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


# cgm_merged.to_csv("dm_proj1_data/intermediate_files/cgm_merged_2.csv")


for i in range(len(cgm_merged)):
    x, y = str(i) + "_x", str(i) + "_y"
    plt.plot(cgm_merged[x], cgm_merged[y])
    plt.title("Patient: 1 Meal #{}".format(x.split("_")[0]))
    plt.xticks(rotation=90)
    plt.savefig("Meal_{}.png".format(x.split("_")[0]))
    plt.close()

num = 2
for i in range(len(cgm_merged.columns) // 2):
    cgm_merged.insert(num, '{}_out_of_range'.format(i), [70 - j if j < 70 else j - 180 if j > 180 else 0 for j in cgm_merged['{}_y'.format(i)]], True)
    num += 3
cgm_merged


posCount, negCount = 0, 0
outOfRange = []
for i in range(len(cgm_merged.columns) // 3):
    for j in cgm_merged['{}_out_of_range'.format(i)]:
        if j > 0:
            posCount += 1
        if j < 0:
            negCount += 1
#     outOfRange.append([posCount,negCount])
    outOfRange.append(posCount)
    posCount, negCount = 0, 0

timeRange = ['{1} - {0}'.format(cgm_merged['{}_x'.format(i)][0], cgm_merged['{}_x'.format(i)][-1]) for i in range(len(cgm_merged.columns) // 3)]
feature_matrix = pd.DataFrame({'timeRange': timeRange, 'outOfRange': outOfRange})

import numpy as np
from scipy.stats import norm, mstats


def mk_test(x, alpha=0.3):
    """   
    Input:
        x:   a vector of data
        alpha: significance level (0.05 default)

    Output:
        trend: tells the trend (increasing, decreasing or no trend)
        h: True (if trend is present) or False (if trend is absence)
        p: p value of the significance test
        z: normalized test statistics 

    Examples
    --------
      >>> x = np.random.rand(100)
      >>> trend,h,p,z = mk_test(x,0.05) 
    """
    n = len(x)

    # calculate S
    s = 0
    for k in range(n - 1):
        for j in range(k + 1, n):
            s += np.sign(x[j] - x[k])

    # calculate the unique data
    unique_x = np.unique(x)
    g = len(unique_x)

    # calculate the var(s)
    if n == g:  # there is no tie
        var_s = (n * (n - 1) * (2 * n + 5)) / 18
    else:  # there are some ties in data
        tp = np.zeros(unique_x.shape)
        for i in range(len(unique_x)):
            tp[i] = np.sum(x == unique_x[i])
        var_s = (n * (n - 1) * (2 * n + 5) - np.sum(tp * (tp - 1) * (2 * tp + 5))) / 18

    if s > 0:
        z = (s - 1) / np.sqrt(var_s)
    elif s < 0:
        z = (s + 1) / np.sqrt(var_s)
    else:  # s == 0:
        z = 0

    # calculate the p_value
    p = 2 * (1 - norm.cdf(abs(z)))  # two tail test
    h = abs(z) > norm.ppf(1 - alpha / 2)

    if (z < 0) and h:
        trend = 'decreasing'
    elif (z > 0) and h:
        trend = 'increasing'
    else:
        trend = 'no trend'

    return trend, h, p, z


# Range length depends on cgm_merged length, if adding a new column please change here
for i in range(len(cgm_merged.columns) // 3):
    col = str(i) + '_y'
    print("Meal #{} :{}".format(i, mk_test(cgm_merged[col].values.ravel())))

# Insert mk_testt_trend
feature_matrix.insert(2, "mk_test_trend", [mk_test(cgm_merged[str(i) + '_y'].values.ravel())[0] for i in range(len(cgm_merged.columns) // 3)])

# Insert mk_test_p_value: p value of the significance test
feature_matrix.insert(3, "mk_test_p_value", [mk_test(cgm_merged[str(i) + '_y'].values.ravel())[2] for i in range(len(cgm_merged.columns) // 3)])

# Insert mk_test_z_value: normalized test statistics
feature_matrix.insert(4, "mk_test_z_value", [mk_test(cgm_merged[str(i) + '_y'].values.ravel())[3] for i in range(len(cgm_merged.columns) // 3)])

feature_matrix.head()
