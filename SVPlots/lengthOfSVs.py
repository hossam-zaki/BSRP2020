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
#
labels = parser.load_obj('labels')


for gene in labels:
    print(gene)
    try:
        geneDF, chromsome = parser.getGeneDF(gene)

        geneData = []

        for index, row in geneDF.iterrows():
            if (int(row['seqnames']) == chromsome) and (int(row['altchr']) == chromsome):
                geneData.append(int(row['altpos']) - int(row['start']))

        data.append(geneData)
    except:
        continue

ind = np.arange(start=0, stop=len(data)*1.5, step=3)
width = 1.25  # the width of the bars
fig2, ax2 = plt.subplots(figsize=(20, 15))
means = []
stdevs = []
print(stats.f_oneway(*data))
for dataset in data:
    means.append(mean(dataset))
for dataset in data:
    try:
        stdevs.append(stdev(dataset))
    except:
        stdevs.append(0)
rects1 = ax2.bar(ind, means, width,
                 label='WT Gene', yerr=stdevs, zorder=0)
for i in range(0, len(ind)):
    # distribute scatter randomly across whole width of bar
    ax2.scatter((ind[i]) + np.random.rand(len(data[i])) *
                width, data[i], color='black', zorder=10, s=10)
ax2.set_ylabel('Number of SVs')
ax2.set_title('Mean number of SVs in WT vs Mutant Genes')
ax2.set_xticks(ind)
ax2.set_xticklabels(labels)
plt.xticks(rotation=90)
plt.scatter([], [], marker=r'$\ast$', label="p < .05",
            color='black', linestyle='None')
plt.scatter([], [], marker=r'$\ast\ast$', label="p < .01",
            color='black', linestyle='None', s=600)
plt.scatter([], [], marker=r'$\ast\ast\ast$',
            label="p < .005", color='black', linestyle='None', s=750)
ax2.legend()
autolabel(rects2, pValues)
plt.show()


plt.show()
