import json
import urllib

import matplotlib.pyplot as plt
import pandas as pd

import helperFiles.buildPlot as plotBuilder

startAndEndJson = json.load(urllib.request.urlopen(
    "https://dcc.icgc.org/api/v1/genes/ENSG00000182185"))  # lift this
start = plotBuilder.lift(startAndEndJson['start'])
end = plotBuilder.lift(startAndEndJson['end'])

df = pd.read_csv('structural_variance_plot/merged_1.6.1.csv')

df = df[(df['seqnames'] == 14) & (df['altchr'] == 14)]

df = df[(df['start'].between(start, end, inclusive=True)) |
        (df['altpos'].between(start, end, inclusive=True))]

unique_ids = df['donor_unique_id'].unique()

print(len(unique_ids))

forPlot = []
for uniID in unique_ids:
    place = df[(df['donor_unique_id'] == uniID)]
    forPlot.append(len(place.index))

fig1, ax1 = plt.subplots()
ax1.set_title('Basic Plot')
ax1.boxplot(forPlot)

plt.show()
