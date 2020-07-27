import json
import urllib
from statistics import mean, median, stdev

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Polygon
from scipy import stats

import dataParser as parser
import helperFiles.buildPlot as plotBuilder

df = pd.read_csv('../merged_1.6.1.csv')

data = []
labels = []

labels = parser.load_obj('labels')
data = parser.load_obj('lengthData')

print(stats.f_oneway(*data))

for i in range(1, len(labels)):
    label = labels[i]
    print(
        f"p-value for RAD51B and {label}: {stats.ttest_ind(data[20], data[i])}")

ind = np.arange(start=0, stop=len(data) * 2, step=2)
width = 1.5  # the width of the bars
fig2, ax2 = plt.subplots(figsize=(20, 15))
means = []
stdevs = []
for dataset in data:
    if len(dataset) == 0:
        means.append(0)
        continue
    means.append(mean(dataset))
for dataset in data:
    try:
        stdevs.append(stdev(dataset))
    except:
        stdevs.append(0)
print(len(ind))
print(len(means))
rects1 = ax2.bar(ind, means, width, yerr=stdevs, zorder=0)
for i in range(0, len(ind)):
    # distribute scatter randomly across whole width of bar
    ax2.scatter((ind[i]) + np.random.rand(len(data[i])) *
                width - width/2, data[i], color='black', zorder=10, s=10)
ax2.set_ylabel('Length in bp')
ax2.set_xlabel('Gene')
ax2.set_title('Mean Length of Intrachromosomal Rearrangements')
ax2.set_xticks(ind)
ax2.set_xticklabels(labels)
plt.xticks(rotation=90)
plt.scatter([], [], marker='.', label="Length of an individual SV",
            color='black', linestyle='None')
ax2.legend()

plt.show()
forBarPlot = []
for gene in labels:
    print(gene)
    try:
        geneDF, chromsome = parser.getGeneDF(gene)

        geneData = 0

        for index, row in geneDF.iterrows():
            if (int(row['seqnames'])) != (int(row['altchr'])):
                geneData += 1
        forBarPlot.append(geneData)
    except:
        forBarPlot.append(0)
        continue

print(forBarPlot)

fig1, ax = plt.subplots(figsize=(10, 6))
ax.set_title('Basic Plot')
x = np.arange(len(labels))
rects = ax.bar(x, forBarPlot)
ax.set_ylabel('Number of Rearragements')
ax.set_xlabel('Gene')
ax.set_title('Number of Interchromosomal Rearrangements per gene')
ax.set_xticks(x)
plt.xticks(rotation=90)
ax.set_xticklabels(labels)

plt.show()
