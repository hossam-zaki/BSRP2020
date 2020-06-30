import json
import re
import subprocess
import sys
import urllib

import ensembl_rest
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import pandas as pd

import helperFiles.buildplot as plotBuilder
import helperFiles.getUniprotRanges as proteinRanges


def getProteinCodingRanges(start, end):
    result = subprocess.run(
        ["Rscript", "protein_domain_plot/getDomains.R", "--start", start, "--end", end], stdout=subprocess.PIPE, text=True)
    str = result.stdout.split('\n')
    toReturn = set()
    for sub in str:
        match = re.match(".*\[\d*\]\s*\d*.(\d*)-(\d*).*", sub)
        if match:
            toReturn.add((match[1], match[2]))
    return toReturn


def buildRanges():
    transcriptAndRange = {}
    ranges = proteinRanges.getRange()
    for feat in ranges:
        print(ranges[feat])
        result = getProteinCodingRanges(ranges[feat][0], ranges[feat][1])
        transcriptAndRange[feat] = []
        for item in result:
            transcriptAndRange[feat].append(item)
    return transcriptAndRange


def buildPlotWithProtein():
    print('building plot...')
    plotBuilder.buildPlot()
    print("loading json...")
    startAndEndJson = json.load(urllib.request.urlopen(
        "https://dcc.icgc.org/api/v1/genes/ENSG00000182185"))  # lift this
    start = plotBuilder.lift(startAndEndJson['start'])
    end = plotBuilder.lift(startAndEndJson['end'])
    print(start)
    ranges = buildRanges()
    cmap = plotBuilder.get_cmap(len(ranges)+1)
    counter = 1
    legend = []
    ax = plt.gca()
    for feat in ranges:
        for transcript in ranges[feat]:
            print(feat)
            print((int(transcript[0])-start),
                  (int(transcript[1])-start))
            rect = patches.Rectangle(
                (int(transcript[0])-start, -counter/500), height=1/500,
                width=((int(transcript[1])-start) - (int(transcript[0])-start)), fill=True, color=cmap(counter-1))
            ax.add_patch(rect)

        legend.append(patches.Patch(color=cmap(counter-1), label=feat))
        counter += 1
    axes = plt.gca()
    axes.set_xlim([0, end-start])
    axes.set_ylim([-counter/500, .05])
    plt.legend(handles=legend,)
    plt.show()


if __name__ == "__main__":
    buildPlotWithProtein()
