import json
import os
import re
import urllib

import pandas as pd

import dataParser

df = pd.read_csv('../merged_1.6.1.csv')
patients = df['donor_unique_id'].unique()
rangeDict = dataParser.buildRangeDict()
if(os.path.isfile("obj/completedDonors0.pkl")):
    completed = dataParser.load_obj("completedDonors0")
else:
    completed = set()
size = int(len(patients)/4)
for patient in patients[0:size+1]:
    if(os.path.isfile("obj/GenewithDonorsWithSVsInGene0.pkl")):
        results = dataParser.load_obj("GenewithDonorsWithSVsInGene0")
    else:
        results = {}
    match = re.match('.*::(.+)', patient)
    donorid = dataParser.getKeyword(match[1])
    if donorid in completed:
        continue
    print(donorid)
    mutsinDonor = dataParser.mutationsInDonorCount(donorid)
    to_nearest_hunderd = 101 - (mutsinDonor % 100)
    try:
        for i in range(0, mutsinDonor+to_nearest_hunderd, 100):
            print(i)
            mutationResponse = json.load(urllib.request.urlopen(
                f"https://dcc.icgc.org/api/v1/donors/{donorid}/mutations?filters=%7B%7D&from={i}&size=100&sort=affectedDonorCountFiltered&order=desc"))
            for j in mutationResponse['hits']:
                if j['type'] == 'single base substitution':
                    chromosome = j['chromosome']
                    mutRange = range(j['start'], j['end'])
                    for gene in rangeDict:
                        if str(rangeDict[gene][0]) == str(chromosome):
                            if dataParser.range_subset(mutRange, rangeDict[gene][1]):
                                if gene not in results:
                                    results[gene] = set()
                                results[gene].add(donorid)
        completed.add(donorid)
        dataParser.save_obj(completed, "completedDonors0")
        dataParser.save_obj(results, "GenewithDonorsWithSVsInGene0")
    except:
        continue
