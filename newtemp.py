from datetime import datetime, timedelta
import pandas as pd
import pandas
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm


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

# for i in range(1, numFiles + 1):
#     cgm_merged.to_csv("dm_proj1_data/intermediate_files/cgm_merged_{}.csv".format(i))


# cgmM = cgm_merged[0]
# for cgmM in cgm_merged:
#     for i in range(len(cgmM)):
#         x, y = str(i) + "_x", str(i) + "_y"
#         plt.plot(cgmM[x], cgmM[y])
#         plt.title("Patient: 1 Meal #{}".format(x.split("_")[0]))
#         plt.xticks(rotation=90)
#         plt.savefig("Meal_{}.png".format(x.split("_")[0]))
#         plt.close()


def addOutofRange(cgmMO):
    num = 2
    for i in range(len(cgmMO.columns) // 2):
        cgmMO.insert(num, '{}_out_of_range'.format(i), [70 - j if j < 70 else j - 180 if j > 180 else 0 for j in cgmMO['{}_y'.format(i)]], True)
        num += 3


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
    outOfRangeList.append(countOutOfRange(cgm_merged[i]))


timeRange = [
    [
        '{1} - {0}'.format(cgm_merged[idx]['{}_x'.format(i)][0], cgm_merged[idx]['{}_x'.format(i)][-1]) for i in range(len(cgm_merged[idx].columns) // 3)
    ] for idx in range(numFiles)
]
feature_matrix = [pd.DataFrame({'timeRange': timeRange[i], 'outOfRange': outOfRangeList[i]}) for i in range(numFiles)]


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

    trend = 'decreasing' if z < 0 and h else 'increasing' if z > 0 and h else 'no trend'

    return trend, h, p, z


# Range length depends on cgm_merged length, if adding a new column please change here
for j in range(numFiles):
    for i in range(len(cgm_merged[j].columns) // 3):
        col = str(i) + '_y'
        print("Meal #{} :{}".format(i, mk_test(cgm_merged[j][col].values.ravel())))

for idx in range(numFiles):
    # Insert mk_testt_trend
    feature_matrix[idx].insert(2, "mk_test_trend", [mk_test(cgm_merged[idx][str(i) + '_y'].values.ravel())[0] for i in range(len(cgm_merged[idx].columns) // 3)])

    # Insert mk_test_p_value: p value of the significance test
    feature_matrix[idx].insert(3, "mk_test_p_value", [mk_test(cgm_merged[idx][str(i) + '_y'].values.ravel())[2] for i in range(len(cgm_merged[idx].columns) // 3)])

    # Insert mk_test_z_value: normalized test statistics
    feature_matrix[idx].insert(4, "mk_test_z_value", [mk_test(cgm_merged[idx][str(i) + '_y'].values.ravel())[3] for i in range(len(cgm_merged[idx].columns) // 3)])


feature_matrix.head()


timeRange = [
    [
        '{1} - {0}'.format(cgm_merged[idx]['{}_x'.format(i)][0], cgm_merged[idx]['{}_x'.format(i)][-1]) for i in range(len(cgm_merged[idx].columns) // 2) if not isNaN(skipCols, idx, i)
    ] for idx in range(numFiles)
]
feature_matrix = [pd.DataFrame({'timeRange': timeRange[i], '#OutOfRange': outOfRangeList[i][0], 'meanOutOfRange':outOfRangeList[i][1]}) for i in range(numFiles)]


for idx in range(numFiles):
    # Insert mk_testt_trend
    feature_matrix[idx].insert(3, "mk_test_trend", [mk_test(cgm_merged[idx][str(i) + '_y'].values.ravel())[0] for i in range(len(cgm_merged[idx].columns) // 2) if not isNaN(skipCols, idx, i)])

    # Insert mk_test_p_value: p value of the significance test
    feature_matrix[idx].insert(4, "mk_test_p_value", [mk_test(cgm_merged[idx][str(i) + '_y'].values.ravel())[2] for i in range(len(cgm_merged[idx].columns) // 2) if not isNaN(skipCols, idx, i)])

    # Insert mk_test_z_value: normalized test statistics
    feature_matrix[idx].insert(5, "mk_test_z_value", [mk_test(cgm_merged[idx][str(i) + '_y'].values.ravel())[3] for i in range(len(cgm_merged[idx].columns) // 2) if not isNaN(skipCols, idx, i)])


for idx in range(numFiles):
    # Insert mk_testt_trend
    print(len(['_x'.format(i) for i in range(len(cgm_merged[idx].columns) // 2) if not isNaN(skipCols, idx, i)]))

    # # Insert mk_test_p_value: p value of the significance test
    # feature_matrix[idx].insert(4, "mk_test_p_value", [mk_test(cgm_merged[idx][str(i) + '_y'].values.ravel())[2] for i in range(len(cgm_merged[idx].columns) // 2) if not isNaN(skipCols, idx, i)])

    # # Insert mk_test_z_value: normalized test statistics
    # feature_matrix[idx].insert(5, "mk_test_z_value", [mk_test(cgm_merged[idx][str(i) + '_y'].values.ravel())[3] for i in range(len(cgm_merged[idx].columns) // 2) if not isNaN(skipCols, idx, i)])

x = StandardScaler().fit_transform(feature_matrix[0][['OutOfRange', 'meanOutOfRange', 'mk_test_trend', 'mk_test_p_value', 'mk_test_z_value', 'variance_FFT', 'MeanRange']])
pca = decomposition.PCA(n_components=7)
pca2 = pca.fit(x)
print(pca2.components_)
print(pca2.explained_variance_ratio_)
# principalComponents = pca.fit_transform(x)
principalComponents = pca2.transform(x)
columns1 = ['OutOfRange', 'meanOutOfRange', 'mk_test_trend', 'mk_test_p_value', 'mk_test_z_value', 'variance_FFT', 'MeanRange']
principalDataframe = pd.DataFrame(data=principalComponents, columns=['OutOfRange', 'meanOutOfRange', 'mk_test_trend', 'mk_test_p_value', 'mk_test_z_value', 'variance_FFT', 'MeanRange'])
#targetDataframe = df[['target']]

#newDataframe = pd.concat([principalDataframe, targetDataframe],axis = 1)
percent_variance = np.round(pca.explained_variance_ratio_ * 100, decimals=2)
columns = ['PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7']
plt.bar(x=range(1, 8), height=percent_variance, tick_label=columns)
plt.ylabel('Percentate of Variance Explained')
plt.xlabel('Principal Component')
plt.title('PCA Scree Plot')
plt.show()