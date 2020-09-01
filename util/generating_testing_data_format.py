"""
Only for preprocessing data.
"""
import pandas as pd
from extract_features_ori import extract_features_pseudoknot_free
from rna_toolkit import bp_to_dp
import os
import csv
import json
import RNA


def RNAfold_str_generation(seq):
    print('*Roll-Out*')
    fc = RNA.fold_compound(OriginalRNAChainCon)

    (s, mm) = fc.mfe()
    print "%s\n%s (MM: %d)\n" %  (seq, s, -mm)
    RNAfold_str = s
    ##print ("The dot-bracket str is",s)
    #break the str s into pieces, and translate it into BPchain
    ####print('show s',s)
    #new_s = [s[i] for i in range(0, len(s))]
    ####print('new_s',new_s)
    #newBP = DP_to_OurBP(new_s)
    return RNAfold_str


df_temp = pd.read_csv(os.path.join(path_of_ori_data, file_name), skiprows=[0], delim_whitespace = True, header=None)

ori_seq = df_temp.iloc[:,1].values.tolist()
seq = ''.join(ori_seq)
print('seq',seq)

sec_str_bp = df_temp.iloc[:,4]
sec_str = bp_to_dp(sec_str_bp)
print('Actual_str',sec_str)

RNAfold_str = RNAfold_str_generation(seq)
print('RNAfold_str',RNAfold_str)



