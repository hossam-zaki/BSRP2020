from __future__ import absolute_import, print_function

import json
import urllib

import icgc
import matplotlib.pyplot as plt
import pandas as pd

import getRangesFromUpdatedEnsembl


def run():
    numberOfSeqs = json.load(urllib.request.urlopen(
        "http://dcc.icgc.org/api/v1/genes/ENSG00000182185/mutations/count"))
    print(numberOfSeqs)
    startAndEndJson = json.load(urllib.request.urlopen(
        "https://dcc.icgc.org/api/v1/genes/ENSG00000182185"))
    start = startAndEndJson['start']
    end = startAndEndJson['end']

    counter = 0
    to_nearest_hunderd = 101 - (numberOfSeqs % 100)
    placeInGenome = []
    numberOfOccurences = []
    for i in range(0, numberOfSeqs+to_nearest_hunderd, 100):
        try:
            json_file = json.load(urllib.request.urlopen(
                f"https://dcc.icgc.org/api/v1/genes/ENSG00000182185/mutations?from={i}&size=100"))
        except:
            print(f"{i} got messed")
        print(i)

        for hit in json_file['hits']:
            if hit['type'] == 'single base substitution':
                placeInGenome.append(hit['start']-start)
                numberOfOccurences.append(hit['affectedDonorCountTotal'])
    plt.stem(placeInGenome, numberOfOccurences)
    plt.xlim((0, end-start))

    ranges = getRangesFromUpdatedEnsembl.getTranscriptsAndRanges()
    xTicks = []
    labels = []
    for key, value in ranges.items():
        labels.append(key)
        for num in value:
            if len(xTicks) == 0:
                xTicks.append(num)
            elif xTicks[len(xTicks) - 1] != num:
                xTicks.append(num)
    plt.xticks(xTicks, labels)
    plt.title("RAD51B Mutation Map")
    plt.xlabel('Position')
    plt.ylabel('# of Donors')
    plt.show()


if __name__ == '__main__':
    run()
