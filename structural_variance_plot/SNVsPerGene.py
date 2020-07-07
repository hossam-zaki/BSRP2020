import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon

import parseCSV as parser

samples = {
    "Nonhomologous end-joining": ["XRCC6",
                                  "XRCC5",
                                  "PRKDC",
                                  "LIG4",
                                  "XRCC4",
                                  "DCLRE1C",
                                  "NHEJ1"],

    "Microhomology end-joining": ["MRE11",
                                  "RAD50",
                                  "NBN",
                                  "RBBP8",
                                  "ERCC4",
                                  "ERCC1",
                                  "LIG1",
                                  "POLL",
                                  "POLB",
                                  "PARP1",
                                  "LIG3",
                                  "XRCC1"
                                  ],
    "Homologous recombination": ["RAD51",
                                 "RAD51B",
                                 "RAD51D",
                                 "DMC1",
                                 "XRCC2",
                                 "XRCC3",
                                 "RAD52",
                                 "RAD54L",
                                 "RAD54B",
                                 "BRCA1",
                                 "SHFM1",
                                 "RAD50",
                                 "MRE11A",
                                 "NBN",
                                 "RBBP8",
                                 "MUS81",
                                 "EME1",
                                 "EME2",
                                 "GIYD1",
                                 "GIYD2",
                                 "GEN1",
                                 ],
    "Base excision repair": ["UNG",
                             "SMUG1",
                             "MBD4",
                             "TDG",
                             "OGG1",
                             "MUTYH",
                             "NTHL1",
                             "MPG",
                             "NEIL1",
                             "NEIL2",
                             "NEIL3",
                             "APEX1",
                             "APEX2",
                             "LIG3",
                             "XRCC1",
                             "PNKP",
                             "APLF",
                             "PARP1",
                             "PARP2",
                             "PARP3",
                             "MGMT",
                             "ALKBH2",
                             "ALKBH3",
                             ],
    "Repair of DNA-topoisomerase crosslinks": ["TDP1",
                                               "TDP2"
                                               ],
    "Mismatch excision repair": ["MSH2"
                                 "MSH3",
                                 "MSH6",
                                 "MLH1",
                                 "PMS2",
                                 "MSH4",
                                 "MSH5",
                                 "MLH3",
                                 "PMS1",
                                 "PMS2L3"
                                 ],
    "Nucleotide excision repair": ["RAD23B",
                                   "CETN2",
                                   "RAD23A",
                                   "XPA",
                                   "DDB1",
                                   "DDB2",
                                   "RPA1",
                                   "RPA2",
                                   "RPA3",
                                   "TFIIH",
                                   "ERCC3",
                                   "ERCC2",
                                   "GTF2H1",
                                   "GTF2H2",
                                   "GTF2H3",
                                   "GTF2H4",
                                   "GTF2H5",
                                   "CDK7",
                                   "CCNH",
                                   "MNAT1",
                                   "ERCC5",
                                   "ERCC1",
                                   "ERCC4",
                                   "LIG1",
                                   "ERCC8",
                                   "ERCC6",
                                   "UVSSA",
                                   "XAB2",
                                   "MMS19",
                                   ],
    "Fanconi anemia": ["FANCA"
                       "FANCB"
                       "FANCC"
                       "BRCA2"
                       "FANCD2"
                       "FANCE"
                       "FANCF"
                       "FANCG"
                       "FANCI"
                       "BRIP1"
                       "FANCL"
                       "FANCM"
                       "PALB2"
                       "RAD51C"
                       "BTBD12"
                       "FAAP20"
                       "FAAP24"
                       ],
    "Modulation of nucleotide pools": ["NUDT1",
                                       "DUT",
                                       "RRM2B",
                                       ],
    "DNA polymerases": ["POLB",
                        "POLG",
                        "POLD1",
                        "POLE",
                        "PCNA",
                        "REV3L",
                        "MAD2L2",
                        "REV1L",
                        "POLH",
                        "POLI",
                        "POLQ",
                        "POLK",
                        "POLL",
                        "POLM",
                        "POLN",
                        ],
    "Editing and processing nucleases": ["FEN1",
                                         "FAN1",
                                         "TREX1",
                                         "TREX2",
                                         "EXO1",
                                         "APTX",
                                         "SPO11",
                                         "ENDOV",
                                         ],
    "Ubiquitination and modification": ["UBE2A"
                                        "UBE2B"
                                        "RAD18"
                                        "SHPRH"
                                        "HLTF"
                                        "RNF168"
                                        "SPRTN"
                                        "RNF8"
                                        "RNF4"
                                        "UBE2V2"
                                        "UBE2N"
                                        ],
    "Chromatin Structure and Modification": ["H2AFX",
                                             "CHAF1A",
                                             "SETMAR",
                                             ],
    "Other conserved DNA damage response genes": ["ATR",
                                                  "ATRIP",
                                                  "MDC1",
                                                  "RAD1",
                                                  "RAD9A",
                                                  "HUS1",
                                                  "RAD17",
                                                  "CHEK1",
                                                  "CHEK2",
                                                  "TP53",
                                                  "TP53BP1",
                                                  "RIF1",
                                                  "TOPBP1",
                                                  "CLK2",
                                                  "PER1",
                                                  ]
}

# labels = [*samples]
labels = []
forplot = []
counter = 0
for i in samples:
    for gene in samples[i]:
        print(gene)
        labels.append(gene)
        toAdd = parser.getNumOfSVs(gene)
        if(toAdd == None):
            continue
        forplot.append(toAdd)

# fig, ax1 = plt.subplots()
# ax1.set_title('Basic Plot')
# print(forplot)
# bp = ax1.boxplot(forplot)

# plt.show()


fig, ax1 = plt.subplots(figsize=(10, 2))
fig.canvas.set_window_title('A Boxplot Example')
fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

bp = ax1.boxplot(forplot, notch=0, sym='+', vert=1, whis=1.5, showmeans=True)
plt.setp(bp['boxes'], color='black')
plt.setp(bp['whiskers'], color='black')
plt.setp(bp['fliers'], color='red', marker='+')

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)

# Hide these grid behind plot objects
ax1.set_axisbelow(True)
ax1.set_title(
    'Comparison of IID Bootstrap Resampling Across Five Distributions')
ax1.set_xlabel('Distribution')
ax1.set_ylabel('Value')

# Now fill the boxes with desired colors
box_colors = ['darkkhaki', 'royalblue']
num_boxes = len(forplot)
medians = np.empty(num_boxes)
for i in range(num_boxes):
    box = bp['boxes'][i]
    boxX = []
    boxY = []
    for j in range(5):
        boxX.append(box.get_xdata()[j])
        boxY.append(box.get_ydata()[j])
    box_coords = np.column_stack([boxX, boxY])
    # Alternate between Dark Khaki and Royal Blue
    ax1.add_patch(Polygon(box_coords, facecolor=box_colors[i % 2]))
    # Now draw the median lines back over what we just filled in
    med = bp['medians'][i]
    medianX = []
    medianY = []
    for j in range(2):
        medianX.append(med.get_xdata()[j])
        medianY.append(med.get_ydata()[j])
        ax1.plot(medianX, medianY, 'k')
    medians[i] = medianY[0]
    # Finally, overplot the sample averages, with horizontal alignment
    # in the center of each box
    ax1.plot(np.average(med.get_xdata()), np.average(forplot[i]),
             color='w', marker='*', markeredgecolor='k')

# Set the axes ranges and axes labels
ax1.set_xlim(0.5, num_boxes + 0.5)
top = 12
bottom = -5
ax1.set_ylim(bottom, top)
ax1.set_xticklabels(labels,
                    rotation=45, fontsize=8)

# Due to the Y-axis scale being different across samples, it can be
# hard to compare differences in medians across the samples. Add upper
# X-axis tick labels with the sample medians to aid in comparison
# (just use two decimal places of precision)
pos = np.arange(num_boxes) + 1
upper_labels = [str(np.round(s, 2)) for s in medians]
weights = ['bold', 'semibold']
for tick, label in zip(range(num_boxes), ax1.get_xticklabels()):
    k = tick % 2
    ax1.text(pos[tick], .95, upper_labels[tick],
             transform=ax1.get_xaxis_transform(),
             horizontalalignment='center', size='x-small',
             weight=weights[k], color=box_colors[k])

# Finally, add a basic legend
fig.text(0.80, 0.08, f'1 Random Numbers',
         backgroundcolor=box_colors[0], color='black', weight='roman',
         size='x-small')
fig.text(0.80, 0.045, 'IID Bootstrap Resample',
         backgroundcolor=box_colors[1],
         color='white', weight='roman', size='x-small')
fig.text(0.80, 0.015, '*', color='white', backgroundcolor='silver',
         weight='roman', size='medium')
fig.text(0.815, 0.013, ' Average Value', color='black', weight='roman',
         size='x-small')

plt.show()
