import json
import urllib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Polygon
from scipy import stats

import dataParser as parser
import helperFiles.buildPlot as plotBuilder


def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='center', color="black", fontsize=6, fontweight="bold")


samples = parser.samples
s = set()
labels = []
numDonors = []
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
            numDonors.append(toAdd)

numSVs = []
counter = 0
s = set()
for i in samples:
    for gene in samples[i]:
        print(gene)
        if(gene in s):
            continue
        toAdd = parser.getNumOfSVs(gene)
        if(toAdd == None):
            continue
        else:
            s.add(gene)
            numSVs.append(toAdd)

fig, ax = plt.subplots(figsize=(15, 10))
ax.set_title('Number of Donors with SVs in Genes and Number of SVs in Genes')
x = np.arange(len(labels))
rects = ax.bar(x, numDonors)
rects1 = ax.bar(x, numSVs, bottom=numDonors)
ax.set_ylabel('Number of Donors/SVs')
ax.set_xlabel('Genes')
ax.set_title('Number of Donors with SVs in Genes and Number of SVs in Genes')
ax.set_xticks(x)
plt.xticks(rotation=90)
ax.set_xticklabels(labels)
plt.legend((rects[0], rects1[0]),
           ('Number of Donors With SV in Gene', 'Number of SVs in Gene'))
# autolabel(rects)
# autolabel(rects1)
plt.show()
