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
        toAdd = parser.getNumOfDonors(gene)
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
ax.set_ylabel('Number of Donors')
ax.set_xlabel('Gene')
ax.set_title('Number of Donors with SVs in Genes')
ax.set_xticks(x)
plt.xticks(rotation=90)
ax.set_xticklabels(labels)
ax.legend()

autolabel(rects)
plt.show()


fig.tight_layout()

plt.show()
