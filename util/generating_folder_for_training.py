"""
dataframe
'file_name', 'seq', 'fe_per','ent_3', 'gc_perentage', 'ensemble_diversity', 'expected_accuracy'

for each RNA name, we generate a csv containing all seqs with the same str:
'dp', 'ensemble_diversity', 'ent_3', 'expected_accuracy', 'fe_per','gc_perentage', 'seq'

"""

import pandas as pd
from extract_features_ori import extract_features_pseudoknot_free
from rna_toolkit import bp_to_dp
import os
import csv
import json

path_of_ori_data = '/home/menghan/Dropbox (ASU)/RNA/ENTRNA-master_new_MH/ENTRNA-master/util/Neg_label'

def Put_something_into_csv(something):
    #python2 can use file instaed of open
    with file("seed_structures_40000_pseudoknot_free_feature.csv", "a") as csvfile: 
        writer = csv.writer(csvfile)
        writer.writerow(something)


#df = pd.DataFrame(columns=['file_name', 'seq', 'ensemble_diversity','ent_3', 'expected_accuracy', 'fe_per', 'gc_perentage'])

dfList = []

for file_item in os.listdir(path_of_ori_data):
    if file_item.endswith('.csv'):
        print('file_item_name:', file_item)
        file_item = os.path.join(path_of_ori_data, file_item)
        df_temp = pd.read_csv(file_item, header=None)
        #print(df_temp)
        #looking for begin of each RNA
        for row in range(0,len(df_temp)):
            #print('this row', df_temp.iloc[row,0])
            if df_temp.iloc[row,0].endswith('.ct'):
                print('Begin')
                print(row+4)
                #print(df_temp.iloc[row+4,0])
                if row+4<len(df_temp):
                    dfList = []
                    file_name = df_temp.iloc[row,0]
                    file_name = file_name[:-3]
                    print('str filename:',file_name)
                    sec_str_dp = df_temp.iloc[row+1,0]
                    bg=3
                    print('row+bg',row+bg)
                    while not df_temp.iloc[row+bg,0].endswith('.ct'):
                        seq = df_temp.iloc[row+bg,0]
                        seq = seq.split()
                        seq = seq[0]
                        #print('seq',seq)

                        df_features = extract_features_pseudoknot_free(seq, sec_str_dp)

                        df_features.insert(loc=0, column='dp',value=sec_str_dp)
                        df_features.insert(loc=1, column='seq',value=seq)
                        #print('df_features',df_features)
                        dfList.append(df_features)
                        bg = bg+1

                        if row+bg >= len(df_temp):
                            break
                    print('Output')
                    df = pd.concat(dfList)
                    df.reset_index(drop=True)
                    csv_file_name = "/home/menghan/Dropbox (ASU)/RNA/ENTRNA-master_new_MH/ENTRNA-master/util/Neg_label_results/{}.csv".format(file_name)
                    df.to_csv(csv_file_name)







