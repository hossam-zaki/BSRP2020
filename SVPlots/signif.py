
import json
import pickle
import urllib
from statistics import mean, median

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Polygon
from scipy import stats

import dataParser as parser
import helperFiles.buildPlot as plotBuilder


def autolabel(rects, significant):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for i in range(0, len(rects)):
        rect = rects[i]
        height = rect.get_height()
        if(significant[i] < (.005/len(significant))):
            annotation = '***'
        elif((significant[i] < (.01/len(significant)))):
            annotation = '**'
        elif((significant[i] < (.05/len(significant)))):
            annotation = '*'
        else:
            annotation = ''
        ax2.annotate(annotation,
                     xy=(rect.get_x() + rect.get_width() / 2, height),
                     xytext=(0, 3),  # 3 points vertical offset
                     textcoords="offset points",
                     ha='center', va='center', color="black", fontsize=12, fontweight='bold')


df = pd.read_csv('../merged_1.6.1.csv')

data = []
labels = []


samples = parser.samples
s = set()
inputs = []
counter = 0
for i in samples:
    for gene in samples[i]:
        if(gene in s):
            continue
        print(gene)
        if(not parser.checkValid(gene)):
            continue
        else:
            s.add(gene)
            inputs.append(gene)
            labels.append(gene)
parser.save_obj(labels, 'labels')
quit()
# pValues = []
# for gene in inputs:
#     print(gene)
#     all_donors = df['donor_unique_id'].unique()
#     donorsInGene = parser.getDonors(gene)
#     diffDonors = parser.diff(all_donors, donorsInGene)

#     affected = []
#     notaffected = []

#     for donor in donorsInGene:
#         donorDF = df[(df['donor_unique_id'] == donor)]
#         affected.append(len(donorDF.index))

#     for donor in diffDonors:
#         donorDF = df[(df['donor_unique_id'] == donor)]
#         notaffected.append(len(donorDF.index))
#     if(len(affected) == 0):
#         affected.append(0)
#     if(len(notaffected) == 0):
#         notaffected.append(0)
#     data.append(notaffected)
#     data.append(affected)
#     labels.append(gene)
#     value = stats.ttest_ind(affected, notaffected)[1]
#     pValues.append(value)
#     print(value)

pValues = parser.load_obj('pvalues')
data = parser.load_obj('data')

# the x locations for the groups
ind = np.arange(start=0, stop=len(data)*1.5, step=3)
width = 1.25  # the width of the bars
fig2, ax2 = plt.subplots(figsize=(20, 15))
means = []
print(ind)
for dataset in data:
    means.append(mean(dataset))
rects1 = ax2.bar(ind-width/2, means[0::2], width,
                 label='WT Gene')
rects2 = ax2.bar(ind + width/2, means[1::2], width,
                 label='MUT Gene')
ax2.set_ylabel('Number of SVs')
ax2.set_title('Mean number of SVs in WT vs Mutant Genes')
ax2.set_xticks(ind)
ax2.set_xticklabels(labels)
plt.xticks(rotation=90)
ax2.legend()
autolabel(rects2, pValues)

plt.show()
