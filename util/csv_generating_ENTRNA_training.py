"""
dataframe
'file_name', 'seq', 'fe_per','ent_3', 'gc_perentage', 'ensemble_diversity', 'expected_accuracy'
"""

import pandas as pd
from extract_features_ori import extract_features_pseudoknot_free
from rna_toolkit import bp_to_dp
import os
import csv
import json

path_of_ori_data = '/home/menghan/Downloads/ENTRNA-master_new_MH/ENTRNA-master/util/seed_structures/'

def Put_something_into_csv(something):
    #python2 can use file instaed of open
    with file("seed_structures_40000_pseudoknot_free_feature.csv", "a") as csvfile: 
        writer = csv.writer(csvfile)
        writer.writerow(something)


#df = pd.DataFrame(columns=['file_name', 'seq', 'ensemble_diversity','ent_3', 'expected_accuracy', 'fe_per', 'gc_perentage'])

dfList = []

itera = 1
for file_name in os.listdir(path_of_ori_data):
    print('the file name is:',file_name)
    print('iteration:',itera)
    if file_name.endswith('.ct'):
        print('Begin')
        df_temp = pd.read_csv(os.path.join(path_of_ori_data, file_name), skiprows=[0], delim_whitespace = True, header=None)
        #print('df_temp is:',df_temp)
        ori_seq = df_temp.iloc[:,1].values.tolist()
        #print('ori_seq',ori_seq)
        seq = ''.join(ori_seq)
        #print(seq)
        sec_str_bp = df_temp.iloc[:,4]
        #print('sec_str_bp',sec_str_bp)
        sec_str = bp_to_dp(sec_str_bp)
        #print('sec_str',sec_str)
        df_features = extract_features_pseudoknot_free(seq, sec_str)

        df_features.insert(loc=0, column='file_name', value=file_name)
        df_features.insert(loc=1, column='seq',value=seq)
        #print('df_features',df_features)
        dfList.append(df_features)
        itera = itera+1
        #Put_something_into_csv(df_features)
    else:
        pass


df = pd.concat(dfList)
df.reset_index(drop=True)
#print(df)
df.to_csv('seed_structures_40000_pseudoknot_free_feature.csv')



