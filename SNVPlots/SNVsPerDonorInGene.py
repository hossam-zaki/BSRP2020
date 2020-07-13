import re

import pandas as pd

import dataParser

df = pd.read_csv('../merged_1.6.1.csv')
patients = df['donor_unique_id'].unique()
for patient in patients:
    match = re.match('.*::(.+)', patient)
    print(dataParser.getKeyword(match[1]))
    quit()
