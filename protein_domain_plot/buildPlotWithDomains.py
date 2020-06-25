import json
import re
import subprocess
import sys
import urllib

import ensembl_rest
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd

import buildPlot as plotBuilder

sys.path.insert(0, 'lollipopplotFiles/buildPlot.py ')


def getTranscriptCodingStartAndEnd(transcriptID, start, width):
    result = subprocess.run(
        ["Rscript", "protein_domain_plot/getDomains.R", "--start", start, "--width", width, "--transcript", transcriptID], stdout=subprocess.PIPE, text=True)
    str = result.stdout.split('\n')
    start = []
    bool = False
    print(result.stdout)
    for sub in str:
        print(bool)
        if bool:
            start.append(re.match(".+\s+(\d+).+", sub)[1])
            bool = False
        if '<integer> <logical>' in sub:
            bool = True
            continue
        match = re.match(f".+{transcriptID}\s+(\d+).*", sub)
        if match:
            start.append(match[1])
    if(len(start) == 2):
        return start
    return -1


def buildRanges():
    json = ensembl_rest.symbol_lookup(
        species='homo sapiens',
        symbol='RAD51B',
        params={'expand': True}
    )
    start = json['start']
    end = json['end']
    transcriptAndRange = {}
    for transcript in json['Transcript']:
        tranStart = transcript['start']
        tranEnd = transcript['end']
        tranID = transcript['id']
        if transcript['biotype'] != 'protein_coding':
            continue
        transcriptAndRange[tranID] = {}
        lengthOfExons = 0
        counter = 1
        for exon in transcript['Exon']:
            exStart = exon['start']
            exEnd = exon['end']
            lengthOfExons += exEnd - exStart + 1
        print(lengthOfExons)
        # for i in range(exStart, exEnd + 1):
        #     transcriptAndRange[tranID][i] = [counter]
        #     counter += 1
        #     x = getTranscriptCodingStartAndEnd(
        #         tranID, str(i), "1")
        #     if x == -1:
        #         transcriptAndRange[tranID][i].append(False)
        #     else:
        #         transcriptAndRange[tranID][i].append(True)
    return transcriptAndRange


def buildPlotWithProtein():

    plotBuilder.buildPlot()
    startAndEndJson = json.load(urllib.request.urlopen(
        "https://dcc.icgc.org/api/v1/genes/ENSG00000182185"))  # lift this
    start = plotBuilder.lift(startAndEndJson['start'])
    end = plotBuilder.lift(startAndEndJson['end'])
    ranges = buildRanges()
    cmap = plotBuilder.get_cmap(len(ranges)+1)
    counter = 1
    legend = []
    for transcript in ranges:
        for nuc in ranges[transcript]:
            if ranges[transcript][nuc][1]:
                plt.axhline(-counter/500, linewidth=8, xmin=(
                    nuc/(end-start)), xmax=(nuc/(end-start)), color=cmap(counter-1))
        legend.append(mpatches.Patch(color=cmap(counter-1), label=transcript))
        counter += 1
    axes = plt.gca()
    axes.set_xlim([0, end-start])
    axes.set_ylim([-counter/500, .05])
    plt.legend(handles=legend,)
    plt.show()


if __name__ == "__main__":
    buildRanges()
