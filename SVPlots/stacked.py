import json
import urllib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Polygon
from scipy import stats

import dataParser as parser
import helperFiles.buildPlot as plotBuilder


def autolabel(rects, bottom):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for i in range(0, len(rects)):
        rect = rects[i]
        height = rect.get_height()
        if height == 0:
            continue
        if len(bottom) == 0:
            val = 0
        else:
            val = bottom[i]
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height + val),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='center', color="black", fontsize=6, fontweight="bold", zorder=10,)


samples = parser.samples
s = set()
labels = []
data = []
counter = 0
for i in samples:
    for gene in samples[i]:
        if(gene in s):
            continue
        print(gene)
        toAdd = parser.getNumOfDonorsWithBuckets(gene)
        if(toAdd == None):
            continue
        else:
            s.add(gene)
            labels.append(gene)
            data.append(np.histogram(toAdd)[0])


maxlen = 0
for arr in data:
    if len(arr) > maxlen:
        maxlen = len(arr)

data = np.array(data)

fig, ax = plt.subplots(figsize=(20, 15))
x = np.arange(len(labels))
bottom = []
forlegend = []
labelsForLegend = []
for i in range(0, maxlen):
    forPlot = []
    for j in data:
        try:
            forPlot.append(j[i])
        except:
            forPlot.append(0)
    if len(bottom) == 0:
        rects = ax.bar(x, forPlot, label=f"{i+1} SVs in Samples")
        #autolabel(rects, bottom)
        bottom = np.array(forPlot)
    else:
        rects = ax.bar(
            x, forPlot, label=f"{i+1} SVs in Samples", bottom=bottom)
        #autolabel(rects, bottom)
        bottom += np.array(forPlot)
    forlegend.append(rects[0])
    labelsForLegend.append(f"{i+1} SVs in Samples")

ax.set_title('Number of Donors with SVs in Samples')
ax.set_ylabel('Number of Donors')
ax.set_xlabel('Genes')
ax.set_title('Number of SVs in Samples')
ax.set_xticks(x)
plt.xticks(rotation=90)
ax.set_xticklabels(labels)
fig.legend(forlegend, labelsForLegend)
plt.show()
