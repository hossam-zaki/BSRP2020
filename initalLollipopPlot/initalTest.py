from __future__ import absolute_import, print_function

import json
import urllib

import icgc
import matplotlib.pyplot as plt
import pandas as pd


def run():
    numberOfSeqs = json.load(urllib.request.urlopen(
        "http://dcc.icgc.org/api/v1/genes/ENSG00000182185/mutations/count"))
    print(numberOfSeqs)
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
            print(hit)
            quit()
            if hit['type'] == 'single base substitution':
                placeInGenome.append(hit['start'])
                numberOfOccurences.append(hit['affectedDonorCountTotal'])
    plt.stem(placeInGenome, numberOfOccurences)
    plt.title("RAD51B Mutation Map")
    plt.xlabel('Position')
    plt.ylabel('# of Donors')
    plt.show()


if __name__ == '__main__':
    run()
