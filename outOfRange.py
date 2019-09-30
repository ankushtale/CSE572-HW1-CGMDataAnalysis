from datetime import datetime, timedelta
import pandas
import matplotlib.pyplot as plt

cgm_timestamps = pandas.read_csv("dm_proj1_data/CGMDatenumLunchPat1.csv")
cgm_values = pandas.read_csv("dm_proj1_data/CGMSeriesLunchPat1.csv")


def mdate_to_pdate(mdate):
    return (datetime.fromordinal(int(mdate)) + timedelta(days=mdate % 1) - timedelta(days=366)).strftime("%Y-%m-%d %H:%M:%S")


# first - last
mdate_to_pdate(737225.584155093) - mdate_to_pdate(737225.48)


for col in cgm_timestamps.columns:
    cgm_timestamps[col] = cgm_timestamps[col].apply(lambda x: mdate_to_pdate(x) if pandas.notnull(x) else x)


cgm_timestamps, cgm_values = cgm_timestamps.T, cgm_values.T

cgm_timestamps = cgm_timestamps.set_index([pandas.Index(['cgm_{}'.format(i) for i in range(1, 32)])])
cgm_values = cgm_values.set_index([pandas.Index(['cgm_{}'.format(i) for i in range(1, 32)])])

cgm_merged = pandas.merge(cgm_timestamps, cgm_values, how="inner", left_index=True, right_index=True)

cols = cgm_merged.columns.tolist()
new_cols = []
for i in range(len(cgm_merged)):
    new_cols.append(cols[0 + i])
    new_cols.append(cols[len(cgm_merged.columns) // 2 + i])


cgm_merged = cgm_merged[new_cols]
cgm_merged = cgm_merged.fillna(method='ffill', limit=2)


cgm_merged.to_csv("dm_proj1_data/intermediate_files/cgm_merged_2.csv")


for i in range(len(cgm_merged)):
    x, y = str(i) + "_x", str(i) + "_y"
    plt.plot(cgm_merged[x], cgm_merged[y])
    plt.title("Meal #{}".format(x.split("_")[0]))
    plt.xticks(rotation=90)
    plt.show()


cgm_merged

num = 2
for i in range(len(cgm_merged.columns) // 2):
    cgm_merged.insert(num, '{}_out_of_range'.format(i), [70 - j if j < 70 else j - 180 if j > 180 else 0 for j in cgm_merged['{}_y'.format(i)]], True)
    num += 3
cgm_merged

# i = 0
# while True:
#     del cgm_merged['{}_out_of_range'.format(i)]
#     i += 1


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


pandas.DataFrame({'timeRange': timeRange, 'outOfRange': outOfRange})
