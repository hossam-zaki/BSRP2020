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


# def autolabel(rects, significant):
#     """Attach a text label above each bar in *rects*, displaying its height."""
#     for i in range(0, len(rects)):
#         rect = rects[i]
#         height = rect.get_height()
#         if(significant[i] < (.005/len(significant))):
#             annotation = '***'
#         elif((significant[i] < (.01/len(significant)))):
#             annotation = '**'
#         elif((significant[i] < (.05/len(significant)))):
#             annotation = '*'
#         else:
#             annotation = ''
#         ax2.annotate(annotation,
#                      xy=(rect.get_x(), 70000000),
#                      xytext=(0, 3),  # 3 points vertical offset
#                      textcoords="offset points",
#                      ha='center', va='center', color="red", fontsize=20, fontweight='bold', zorder=20)
def autolabel(rects, xpos='center'):
    """
    Attach a text label above each bar in *rects*, displaying its height.

    *xpos* indicates which side to place the text w.r.t. the center of
    the bar. It can be one of the following {'center', 'right', 'left'}.
    """

    ha = {'center': 'center', 'right': 'left', 'left': 'right'}
    offset = {'center': 0, 'right': 1, 'left': -1}

    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(offset[xpos]*3, 3),  # use 3 points offset
                    textcoords="offset points",  # in both directions
                    ha=ha[xpos], va='bottom')

# for gene in labels:
#     print(gene)
#     try:
#         geneDF, chromsome = parser.getGeneDF(gene)
#         geneData = []
#         for index, row in geneDF.iterrows():
#             if (int(row['seqnames']) == chromsome) and (int(row['altchr']) == chromsome):
#                 geneData.append(int(row['altpos']) - int(row['start']))
#         if len(geneData) == 0:
#             data.append([0])
#             continue
#         data.append(geneData)
#     except:
#         data.append([0])
#     print(data[len(data) - 1])
# parser.save_obj(data, 'lengthDataZoomed')

# print(stats.f_oneway(*data))
# pValues = []
# for i in range(1, len(labels)):
#     if(i == 2):
#         pValues.append(1.0)
#     print(labels[i])
#     print(stats.ttest_ind(data[2], data[i])[1])
#     pValues.append(stats.ttest_ind(data[2], data[i])[1])
# print(len(data))
# print(len(pValues))
# ind = np.arange(start=0, stop=len(data) * 2, step=2)
# width = 1.5  # the width of the bars
# fig2, ax2 = plt.subplots(figsize=(20, 15))
# means = []
# stdevs = []
# for dataset in data:
#     if len(dataset) == 0:
#         means.append(0)
#         continue
#     means.append(mean(dataset))
# for dataset in data:
#     try:
#         stdevs.append(stdev(dataset))
#     except:
#         stdevs.append(0)
# print(len(ind))
# print(len(means))
# rects1 = ax2.bar(ind, means, width, yerr=stdevs, zorder=0)
# for i in range(0, len(ind)):
#     # distribute scatter randomly across whole width of bar
#     ax2.scatter((ind[i]) + np.random.rand(len(data[i])) *
#                 width - width/2, data[i], color='black', zorder=10, s=10)
# ax2.set_ylabel('Length in bp')
# ax2.set_xlabel('Gene')
# ax2.set_title('Mean Length of Intrachromosomal Rearrangements')
# ax2.set_xticks(ind)
# ax2.set_xticklabels(labels)
# plt.xticks(rotation=90)
# autolabel(rects1, pValues)
# plt.scatter([], [], marker='.', label="Length of an individual SV",
#             color='black', linestyle='None')
# plt.scatter([], [], marker=r'$\ast$', label="p < .05",
#             color='red', linestyle='None')
# plt.scatter([], [], marker=r'$\ast\ast$', label="p < .01",
#             color='red', linestyle='None', s=600)
# plt.scatter([], [], marker=r'$\ast\ast\ast$',
#             label="p < .005", color='red', linestyle='None', s=750)
# ax2.legend()

# plt.show()
# quit()


interchromosomal = []
intrachromosomal = []
samples = parser.samples
labels = []
s = set()
counter = 0
# for category in samples:
#     for gene in samples[category]:
#         if(gene in s):
#             continue
#         print(gene)
#         if(not parser.checkValid(gene)):
#             continue
#         else:
#             s.add(gene)
#             labels.append(gene)
labels = parser.load_obj('WTvsMUTsubsetLabels')
for gene in labels:
    print(gene)
    try:
        geneDF, chromsome = parser.getGeneDF(gene)
        geneDataIntra = 0
        geneDataInter = 0
        for index, row in geneDF.iterrows():
            if (int(row['seqnames']) == chromsome) and (int(row['altchr']) == chromsome):
                geneDataIntra += 1
            if (int(row['seqnames'])) != (int(row['altchr'])):
                geneDataInter += 1
        intrachromosomal.append(geneDataIntra)
        interchromosomal.append(geneDataInter)
    except:
        intrachromosomal.append(0)
        interchromosomal.append(geneDataInter)

fig1, ax = plt.subplots(figsize=(10, 6))

ind = np.arange(len(interchromosomal))  # the x locations for the groups
width = 0.35  # the width of the bars
print(interchromosomal)
fig, ax = plt.subplots()
rects1 = ax.bar(ind - width/2, interchromosomal, width,
                label='Number of Interchromosomal Rearrangements')
rects2 = ax.bar(ind + width/2, intrachromosomal, width,
                label='Number of Intrachromosomal Rearrangements')

ax.set_ylabel('Number of Rearrangements', fontsize=16)
ax.set_xlabel('Genes')
ax.set_xticks(ind)
ax.set_xticklabels(labels, fontsize=13)
plt.xticks(rotation=90)
ax.legend(prop={'size': 13})
autolabel(rects1)
autolabel(rects2)

fig1.tight_layout()
plt.show()
