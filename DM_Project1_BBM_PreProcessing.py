import pandas

insulinBasal = pandas.read_csv('dm_proj1_data/InsulinBasalLunchPat1.csv')
insulinBolus = pandas.read_csv('dm_proj1_data/InsulinBolusLunchPat1.csv')
insulinDatenum = pandas.read_csv('dm_proj1_data/InsulinDatenumLunchPat1.csv')

for col in insulinDatenum.columns:
    insulinDatenum[col] = pandas.to_datetime(insulinDatenum[col], unit='s')

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
