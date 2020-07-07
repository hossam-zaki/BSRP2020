import json
import urllib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Polygon

import helperFiles.buildPlot as plotBuilder


def getNumOfSVs(symb):
    df = pd.read_csv('structural_variance_plot/merged_1.6.1.csv')
    try:
        getSymbol = json.load(urllib.request.urlopen(
            f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symb}?content-type=application/json;expand=1"))
        ensembSymb = getSymbol['id']
        startAndEndJson = json.load(urllib.request.urlopen(
            f"https://dcc.icgc.org/api/v1/genes/{ensembSymb}"))  # lift this
        chromosome = int(startAndEndJson['chromosome'])
        if startAndEndJson['assembly_name'] == 'GRCh37':
            start = plotBuilder.lift(startAndEndJson['start'], chromosome)
            end = plotBuilder.lift(startAndEndJson['end'], chromosome)
        else:
            start = startAndEndJson['start']
            end = startAndEndJson['end']
    except:
        print(f"symb got fricked")
    df = df[(df['seqnames'] == chromosome) | (df['altchr'] == chromosome)]

    df = df[(df['start'].between(start, end, inclusive=True)) |
            (df['altpos'].between(start, end, inclusive=True))]

    unique_ids = df['donor_unique_id'].unique()
    forPlot = []
    for uniID in unique_ids:
        place = df[(df['donor_unique_id'] == uniID)]
        forPlot.append(len(place.index))
    return forPlot


# fig, ax1 = plt.subplots()
# ax1.set_title('Basic Plot')
# data = getNumOfSVs('RAD51B')
# bp = ax1.boxplot(data, notch=0, sym='+', vert=1, whis=1.5, showmeans=True)
# plt.setp(bp['boxes'], color='black')
# plt.setp(bp['whiskers'], color='black')
# plt.setp(bp['fliers'], color='red', marker='+')

# ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
#                alpha=0.5)

# # Hide these grid behind plot objects
# ax1.set_axisbelow(True)
# ax1.set_title(
#     'Comparison of IID Bootstrap Resampling Across Five Distributions')
# ax1.set_xlabel('Distribution')
# ax1.set_ylabel('Value')
# # print(bp['means'][0].get_ydata()[0])
# # quit()
# # Now fill the boxes with desired colors
# box_colors = ['darkkhaki', 'royalblue']
# num_boxes = len(data)
# medians = np.empty(num_boxes)
# for i in range(0, 1):
#     box = bp['boxes'][i]
#     boxX = []
#     boxY = []
#     for j in range(5):
#         boxX.append(box.get_xdata()[j])
#         boxY.append(box.get_ydata()[j])
#     box_coords = np.column_stack([boxX, boxY])
#     # Alternate between Dark Khaki and Royal Blue
#     ax1.add_patch(Polygon(box_coords, facecolor=box_colors[i % 2]))
#     # Now draw the median lines back over what we just filled in
#     med = bp['means'][i]
#     print(med)
#     medianX = []
#     medianY = []
#     for j in range(1):
#         medianX.append(med.get_xdata()[j])
#         medianY.append(med.get_ydata()[j])
#         ax1.plot(medianX, medianY, 'k')
#     medians[i] = medianY[0]
#     # Finally, overplot the sample averages, with horizontal alignment
#     # in the center of each box
#     ax1.plot(np.average(med.get_xdata()), np.average(data[i]),
#              color='w', marker='*', markeredgecolor='k')

# # Set the axes ranges and axes labels
# ax1.set_xlim(0.5, num_boxes + 0.5)
# top = 40
# bottom = -5
# ax1.set_ylim(bottom, top)
# random_dists = ['RAD51B']
# ax1.set_xticklabels(np.repeat(random_dists, 2),
#                     rotation=45, fontsize=8)

# # Due to the Y-axis scale being different across samples, it can be
# # hard to compare differences in medians across the samples. Add upper
# # X-axis tick labels with the sample medians to aid in comparison
# # (just use two decimal places of precision)
# pos = np.arange(num_boxes) + 1
# upper_labels = [str(np.round(s, 2)) for s in medians]
# weights = ['bold', 'semibold']
# for tick, label in zip(range(num_boxes), ax1.get_xticklabels()):
#     k = tick % 2
#     ax1.text(pos[tick], .95, upper_labels[tick],
#              transform=ax1.get_xaxis_transform(),
#              horizontalalignment='center', size='x-small',
#              weight=weights[k], color=box_colors[k])
# N = 1
# # Finally, add a basic legend
# fig.text(0.80, 0.08, f'{N} Random Numbers',
#          backgroundcolor=box_colors[0], color='black', weight='roman',
#          size='x-small')
# fig.text(0.80, 0.045, 'IID Bootstrap Resample',
#          backgroundcolor=box_colors[1],
#          color='white', weight='roman', size='x-small')
# fig.text(0.80, 0.015, '*', color='white', backgroundcolor='silver',
#          weight='roman', size='medium')
# fig.text(0.815, 0.013, ' Average Value', color='black', weight='roman',
#          size='x-small')

# plt.show()
