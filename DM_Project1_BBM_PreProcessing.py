import pandas
import datetime
import matplotlib

insulinBasal = pandas.read_csv('dm_proj1_data/InsulinBasalLunchPat1.csv')
insulinBolus = pandas.read_csv('dm_proj1_data/InsulinBolusLunchPat1.csv')
insulinDatenum = pandas.read_csv('dm_proj1_data/InsulinDatenumLunchPat1.csv')


for col in insulinDatenum.columns:
    insulinDatenum[col] = insulinDatenum[col].apply(lambda x: (datetime.datetime.fromordinal(int(x)) + datetime.timedelta(days=x % 1) - datetime.timedelta(days=366)).strftime('%Y-%m-%d %H:%M:%S') if pandas.notnull(x) else x)


insulinBasal = insulinBasal.transpose()
insulinBolus = insulinBolus.transpose()
insulinDatenum = insulinDatenum.transpose()

insulinBasal = insulinBasal.set_index([pandas.Index(['cgm_{}'.format(i) for i in range(1, 40)])])
insulinBolus = insulinBolus.set_index([pandas.Index(['cgm_{}'.format(i) for i in range(1, 40)])])
insulinDatenum = insulinDatenum.set_index([pandas.Index(['cgm_{}'.format(i) for i in range(1, 40)])])

firstMerge = pandas.merge(insulinDatenum, insulinBasal, how="inner", left_index=True, right_index=True)

secondMerge = pandas.merge(firstMerge, insulinBolus, how="inner", left_index=True, right_index=True)

newCols = []
tempCols = secondMerge.columns.tolist()
for i in range(len(tempCols) // 3):
    newCols.append(tempCols[i])
    newCols.append(tempCols[len(insulinBasal.columns.tolist()) + i])
    newCols.append(tempCols[2 * len(insulinBasal.columns.tolist()) + i])

finalDF = secondMerge[newCols]

finalDF.to_csv('dm_proj1_data/intermediate_files/basal_bolus_merged.csv')
