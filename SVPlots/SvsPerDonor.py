
import json
import urllib.request

import matplotlib.pyplot as plt

numberOfDonors = json.load(urllib.request.urlopen(
    "https://dcc.icgc.org/api/v1/genes/ENSG00000182185/donors/count"))

to_nearest_hunderd = 101 - (numberOfDonors % 100)
donorDict = {}
forScatter = []
y = []
#
counter = 0
for i in range(0, numberOfDonors+to_nearest_hunderd, 100):
    try:
        json_file = json.load(urllib.request.urlopen(
            f"https://dcc.icgc.org/api/v1/genes/ENSG00000182185/donors?filters=%7B%7D&from={i}&size=100"))
    except:
        print(f"{i} got messed")
    print(i)

    for hit in json_file['hits']:
        donID = hit['id']
        numOfMuts = json.load(urllib.request.urlopen(
            f"https://dcc.icgc.org/api/v1/donors/{donID}/genes/ENSG00000182185/mutations/counts"))
        forScatter.append(numOfMuts[donID]["ENSG00000182185"])
        y.append(1)

plt.scatter(y, forScatter)
plt.show()
