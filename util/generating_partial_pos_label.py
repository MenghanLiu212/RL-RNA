"""
This code is to delete all the rows that does not fit file_name in ned_label folder.

"""
import pandas as pd
from extract_features_ori import extract_features_pseudoknot_free
from rna_toolkit import bp_to_dp
import os
import csv
import json

df_temp = pd.read_csv('seed_structures_40000_pseudoknot_free_feature_new_partial.csv')

name_list = df_temp['file_name']

#path_of_ori_data = '/home/menghan/Dropbox (ASU)/RNA/ENTRNA-master_new_MH/ENTRNA-master/util/Neg_label_results'
path_of_ori_data = '/home/menghan/Dropbox (ASU)/RNA/ENTRNA-master_new_MH/ENTRNA-master/util/Neg_label_results_2000'

count_drop=0
drop_list=[]
for i in range(0,len(name_list)):
    count=0
    for file_item in os.listdir(path_of_ori_data):
        if file_item.endswith('.csv'):
            file_name = str(file_item)[:-4]
            if name_list[i] == file_name:
                print('FLAG')
                print('file_name',file_name)
                print('name_list[i]', name_list[i])
                break
            else:
                count=count+1
    if count==2000:
        #df_temp.drop([i])
        print('drop:', name_list[i], 'num', i)
        drop_list.append(i)
        count_drop=count_drop+1




df_temp.drop(df_temp.index[drop_list], inplace=True)
#df_temp.to_csv('seed_structures_40000_pseudoknot_free_feature_partial.csv')
df_temp.to_csv('seed_structures_2000_pseudoknot_free_feature_partial.csv')
print('drop_list', drop_list)
print('count_drop', count_drop)

