from __future__ import absolute_import, print_function

import json
import urllib

import icgc
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
from pyliftover import LiftOver

import getRangesFromUpdatedEnsembl

lo = LiftOver('hg19', 'hg38')


def lift(coord):
    return lo.convert_coordinate('chr14', coord)[0][1]


def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)


def run():
    numberOfSeqs = json.load(urllib.request.urlopen(
        "http://dcc.icgc.org/api/v1/genes/ENSG00000182185/mutations/count"))
    startAndEndJson = json.load(urllib.request.urlopen(
        "https://dcc.icgc.org/api/v1/genes/ENSG00000182185"))  # lift this
    start = lift(startAndEndJson['start'])
    end = lift(startAndEndJson['end'])
    counter = 0
    to_nearest_hunderd = 101 - (numberOfSeqs % 100)
    placeInGenome = []
    numberOfOccurences = []
    #
    for i in range(0, numberOfSeqs+to_nearest_hunderd, 100):
        try:
            json_file = json.load(urllib.request.urlopen(
                f"https://dcc.icgc.org/api/v1/genes/ENSG00000182185/mutations?from={i}&size=100"))
        except:
            print(f"{i} got messed")
        print(i)

        for hit in json_file['hits']:
            if hit['type'] == 'single base substitution':
                placeInGenome.append(lift(hit['start'])-start)  # lift this
                numberOfOccurences.append(
                    hit['affectedDonorCountTotal']/hit['testedDonorCount'])
    plt.stem(placeInGenome, numberOfOccurences)

    ranges = getRangesFromUpdatedEnsembl.getTranscriptsAndRanges()
    cmap = get_cmap(len(ranges)+1)
    counter = 1
    legend = []
    for transcript in ranges:
        for exon in ranges[transcript]:
            print(exon)
            plt.axhline(-counter/500, linewidth=8, xmin=(
                exon[0]/(end-start)), xmax=(exon[1]/(end-start)), color=cmap(counter-1))
        legend.append(mpatches.Patch(color=cmap(counter-1), label=transcript))
        counter += 1
        # line.set_label('yee')
        # for num in value:
        #     if len(xTicks) == 0:
        #         xTicks.append(num)
        #     elif xTicks[len(xTicks) - 1] != num:
        #         xTicks.append(num)
    axes = plt.gca()

    axes.set_ylim([-counter/500, .05])
    plt.legend(handles=legend,)
    plt.title("RAD51B Mutation Map")
    plt.xlabel('Position')
    plt.ylabel('Frequency of Mutation')
    plt.show()


if __name__ == '__main__':
    run()
