
import pandas as pd
from extract_features_ori import extract_features_pseudoknot_free
from rna_toolkit import bp_to_dp
import os
import csv
import json

df_temp = pd.read_csv('seed_structures_40000_pseudoknot_free_feature.csv')

temp = df_temp['file_name']

print(temp)


for i in range(0, len(temp)):
    temp_row = temp[i]
    print(temp_row)
    temp_row = temp_row[:-3]
    temp[i]=temp_row

df_temp.to_csv('seed_structures_40000_pseudoknot_free_feature_new.csv')

