import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon

import dataParser as parser

samples = parser.samples

# labels = [*samples]
labels = []
forplot = []
counter = 0
for i in samples:
    for gene in samples[i]:
        print(gene)
        toAdd = parser.getNumOfSVs(gene)
        if(toAdd == None):
            continue
        else:
            labels.append(gene)
            forplot.append(toAdd)


def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')


fig, ax = plt.subplots()
ax.set_title('Basic Plot')
print(forplot)
x = np.arange(len(labels))
rects = ax.bar(x, forplot)
ax.set_ylabel('Number of Total SVs')
ax.set_title('Number of SV in Genes')
ax.set_xticks(x)
plt.xticks(rotation=90)
ax.set_xticklabels(labels)
ax.legend()

autolabel(rects)
plt.show()


fig.tight_layout()

plt.show()

# fig, ax1 = plt.subplots(figsize=(10, 2))
# fig.canvas.set_window_title('A Boxplot Example')
# fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

# bp = ax1.boxplot(forplot, notch=0, sym='+', vert=1, whis=1.5, showmeans=True)
# plt.setp(bp['boxes'], color='black')
# plt.setp(bp['whiskers'], color='black')
# plt.setp(bp['fliers'], color='red', marker='+')

# # Add a horizontal grid to the plot, but make it very light in color
# # so we can use it for reading data values but not be distracting
# ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
#                alpha=0.5)

# # Hide these grid behind plot objects
# ax1.set_axisbelow(True)
# ax1.set_title(
#     'Comparison of IID Bootstrap Resampling Across Five Distributions')
# ax1.set_xlabel('Distribution')
# ax1.set_ylabel('Value')

# # Now fill the boxes with desired colors
# box_colors = ['darkkhaki', 'royalblue']
# num_boxes = len(forplot)
# medians = np.empty(num_boxes)
# for i in range(num_boxes):
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
#     med = bp['medians'][i]
#     medianX = []
#     medianY = []
#     for j in range(2):
#         medianX.append(med.get_xdata()[j])
#         medianY.append(med.get_ydata()[j])
#         ax1.plot(medianX, medianY, 'k')
#     medians[i] = medianY[0]
#     # Finally, overplot the sample averages, with horizontal alignment
#     # in the center of each box
#     ax1.plot(np.average(med.get_xdata()), np.average(forplot[i]),
#              color='w', marker='*', markeredgecolor='k')

# # Set the axes ranges and axes labels
# ax1.set_xlim(0.5, num_boxes + 0.5)
# top = 12
# bottom = -5
# ax1.set_ylim(bottom, top)
# ax1.set_xticklabels(labels,
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

# # Finally, add a basic legend
# fig.text(0.80, 0.08, f'1 Random Numbers',
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
