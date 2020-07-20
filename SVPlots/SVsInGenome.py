import json
import urllib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Polygon

import dataParser as parser
import helperFiles.buildPlot as plotBuilder

df = pd.read_csv('../merged_1.6.1.csv')

data = []

# for label in parser.samples:
#     for gene in parser.samples[label]:
donors = parser.getDonors('DDB2')
for donor in donors:
    donorDF = df[(df['donor_unique_id'] == donor)]
    data.append(len(donorDF.index))

print(data)
data.sort()

n, bins, patches = plt.hist(x=data, bins='auto', color='#0504aa',
                            alpha=0.7, rwidth=0.85)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('My Very Own Histogram')
plt.text(23, 45, r'$\mu=15, b=3$')
maxfreq = n.max()
# Set a clean upper y-axis limit.
plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)

plt.show()
