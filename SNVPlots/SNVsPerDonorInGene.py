import re

import pandas as pd

import dataParser

df = pd.read_csv('../merged_1.6.1.csv')
patients = df['donor_unique_id'].unique()
for patient in patients:
    match = re.match('.*::(.+)', patient)
    donorid = dataParser.getKeyword(match[1])
    mutsinDonor = dataParser.mutationsInDonorCount(donorid)
    to_nearest_hunderd = 101 - (mutsinDonor % 100)
    for i in range(0, mutsinDonor+to_nearest_hunderd, 100):
