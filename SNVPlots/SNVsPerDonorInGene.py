import json
import re
import urllib

import pandas as pd

import dataParser

df = pd.read_csv('../merged_1.6.1.csv')
patients = df['donor_unique_id'].unique()
rangeDict = dataParser.buildRangeDict()
results = {}
for patient in patients:
    match = re.match('.*::(.+)', patient)
    donorid = dataParser.getKeyword(match[1])
    if donorid == None:
        continue
    print(donorid)
    mutsinDonor = dataParser.mutationsInDonorCount(donorid)
    to_nearest_hunderd = 101 - (mutsinDonor % 100)
    for i in range(0, mutsinDonor+to_nearest_hunderd, 100):
        print(i)
        mutationResponse = json.load(urllib.request.urlopen(
            f"https://dcc.icgc.org/api/v1/donors/{donorid}/mutations?filters=%7B%7D&from={i}&size=100&sort=affectedDonorCountFiltered&order=desc"))
        for j in mutationResponse['hits']:
            if j['type'] == 'single base substitution':
                chromosome = j['chromosome']
                mutRange = range(j['start'], j['end'])
                for gene in rangeDict:
                    if rangeDict[gene][0] == chromosome:
                        if dataParser.range_subset(mutRange, rangeDict[gene][1]):
                            if gene not in results:
                                results[gene] = [donorid]
                            else:
                                results[gene].append(donorid)
dataParser.save_obj(results, "GenewithDonorsWithSVsInGene")
