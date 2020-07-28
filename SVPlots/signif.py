
import json
import pickle
import urllib
from statistics import mean, median, stdev

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Polygon
from scipy import stats

import dataParser as parser
import helperFiles.buildPlot as plotBuilder

#plt.rcParams.update({'font.size': 12})

# plt.xticks(rotation=45)


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
                     xy=(rect.get_x(), 2000),
                     xytext=(0, 3),  # 3 points vertical offset
                     textcoords="offset points",
                     ha='center', va='center', color="red", fontsize=20, fontweight='bold', zorder=20)


df = pd.read_csv('../merged_1.6.1.csv')

data = []
labels = []


samples = parser.samples
s = set()
inputs = []
counter = 0
for gene in samples["Homologous recombination"]:
    if(gene in s):
        continue
    print(gene)
    if(not parser.checkValid(gene)):
        continue
    else:
        s.add(gene)
        inputs.append(gene)
        labels.append(gene)
pValues = []
for gene in inputs:
    print(gene)
    all_donors = df['donor_unique_id'].unique()
    donorsInGene = parser.getDonors(gene)
    diffDonors = parser.diff(all_donors, donorsInGene)

    affected = []
    notaffected = []

    for donor in donorsInGene:
        donorDF = df[(df['donor_unique_id'] == donor)]
        affected.append(len(donorDF.index))

    for donor in diffDonors:
        donorDF = df[(df['donor_unique_id'] == donor)]
        notaffected.append(len(donorDF.index))
    if(len(affected) == 0):
        affected.append(0)
    if(len(notaffected) == 0):
        notaffected.append(0)
    data.append(notaffected)
    data.append(affected)
    value = stats.ttest_ind(affected, notaffected)[1]
    pValues.append(value)
    print(value)

parser.save_obj(pValues, 'WTvsmUTsubsetpValues')
parser.save_obj(data, 'WTvsMUTsubsetData')
parser.save_obj(labels, 'WTvsMUTsubsetLabels')
# pValues = parser.load_obj('WTvsMUTpvalues')
# data = parser.load_obj('WTvsMUTdata')
# labels = parser.load_obj('labels')
# the x locations for the groups
ind = np.arange(start=0, stop=len(data)*1.5, step=3)
width = 1.25  # the width of the bars
fig2, ax2 = plt.subplots(figsize=(20, 15))
means = []
stdevs = []
print(ind)
for dataset in data:
    means.append(mean(dataset))
for dataset in data:
    try:
        stdevs.append(stdev(dataset))
    except:
        stdevs.append(0)
rects1 = ax2.bar(ind-width/2, means[0::2], width,
                 label='WT Gene', yerr=stdevs[0::2], zorder=0)
rects2 = ax2.bar(ind + width/2, means[1::2], width,
                 label='MUT Gene', yerr=stdevs[1::2], zorder=0)
for i in range(0, len(ind)):
    # distribute scatter randomly across whole width of bar

    ax2.scatter((ind[i]-width/2) + np.random.rand(len(data[i + i])) *
                width - width/2, data[i + i], color='blue', edgecolor='black', zorder=10, s=10)
    ax2.scatter((ind[i]+width/2) + np.random.rand(len(data[i + i + 1])) *
                width - width/2, data[i + i + 1], color='orange', edgecolor='black', zorder=10, s=10)
ax2.set_ylabel('Number of SVs')
ax2.set_title('Mean number of SVs in WT vs Mutant Genes')
ax2.set_xticks(ind)
ax2.set_xticklabels(labels)
plt.xticks(rotation=90)
plt.scatter([], [], marker='.', label="Individual Sample with MUT Gene",
            color='orange', edgecolor='black', zorder=10, s=30, linestyle='None')
plt.scatter([], [], marker='.', label="Individual Sample with WT Gene",
            color='blue', edgecolor='black', zorder=10, s=30, linestyle='None')
plt.scatter([], [], marker=r'$\ast$', label="p < .05",
            color='red', linestyle='None')
plt.scatter([], [], marker=r'$\ast\ast$', label="p < .01",
            color='red', linestyle='None', s=600)
plt.scatter([], [], marker=r'$\ast\ast\ast$',
            label="p < .005", color='red', linestyle='None', s=750)
# l = ax2.legend()
# l.set_zorder(50)
autolabel(rects2, pValues)
plt.show()
