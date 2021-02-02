"""
Recalculate all the foldabilities using only-U
This code's path should be inside ENTRNA-master_new_MH.
"""

import RNA
#from ENTRNA import entrna_main
from util.pseudoknot_free import entrna_main
from util.pseudoknot_free_ori import entrna_main_ori
from util.pseudoknot_free import entrna_train_our_ver
#from util.pseudoknot_free import bulid_model
from util.pseudoknot_free import entrna_main_return_all_features
from util.pseudoknot_free_ori import entrna_main_ori
from util.pseudoknot_free_ori import entrna_train_our_ver_ori
from util.pseudoknot_free_ori import entrna_main_return_all_features_ori
import pandas as pd
import numpy
import csv
import json
from pandas import DataFrame
import os
import Main_SeqAttRollOut_RNAfold_only_MFE_4_branch_4_min_Agave_cleared as Our_alg_Main
import time
#*************************************************************

def Train_ENTRNA():
    scaler, clf = entrna_train_our_ver()
    return scaler, clf


def Train_ENTRNA_ori():
    scaler, clf = entrna_train_our_ver_ori()
    return scaler, clf


def Calculate_Foldability(seq1,dp_str,scaler_ori, clf_ori):
    foldability = entrna_main_ori(seq1,dp_str,scaler_ori, clf_ori)
    #foldability = entrna_main(seq1,dp_str,scaler, clf)
    #print('foldability',foldability)
    return foldability


def Calculate_Foldability_NFE(seq, dp_str, scaler, clf):
    #print '*Calculate Foldability NFE'
    #print 'seq:',seq
    #print 'dp_str:', dp_str
    foldability_WOMFE = entrna_main(seq, dp_str, scaler, clf)
    #foldability = entrna_main(seq1,dp_str,scaler, clf)
    #print('foldability',foldability)
    return foldability_WOMFE

def Transfer_Ori_Seq_To_Our_Form(Ori_Seq):
    new_seq = [Ori_Seq[i] for i in range(0,len(Ori_Seq))]
    return new_seq

def Transfer_T_into_U(seq):
    for i in range(len(seq)):
        new_seq = seq
        if new_seq[i] == 'T':
            new_seq[i] = 'U'
    return new_seq

def Read_txt(direct,file_name):
    with open(os.path.join(direct, file_name), "r") as f:
        mylist = f.read().splitlines()
    #print(mylist[0])
    #This_name = mylist[0]
    return mylist

def ModifyingWholeChainByRNAfold(s):
    s = [s[i] for i in range(0, len(s))]
    for i in range(0,len(s)-4):
        if s[i:i+5] == ['(','.','.','.',')']:
            s[i] = '.'
            s[i+4] = '.'
        elif s[i:i+4] == ['(','.','.',')']:
            s[i] = '.'
            s[i + 3] = '.'
        elif s[i:i+3] == ['(','.',')']:
            s[i] = '.'
            s[i + 2] = '.'
    s = ''.join(map(str, s))
    return s

#*************************************************************
path_of_original_csv = '/home/menghan/Dropbox (ASU)/RL-RNA/sanitized_datasets'
path_of_destination = '/home/menghan/Dropbox (ASU)/RNA/Alg_result_analysis/Takcling_ENTRNA_U_to_T_issue/results'
#path_of_txt = '/home/menghan/Dropbox (ASU)/RNA/Alg_result_analysis/Takcling_ENTRNA_U_to_T_issue'
csv_filename = 'rfam_MFENFE.csv'
#txt_filename = 'Matthews_U_name_MFE.txt'

df = pd.read_csv(os.path.join(path_of_original_csv, csv_filename))
#txt_name_list = Read_txt(path_of_txt,txt_filename)

scaler, clf = Train_ENTRNA()
scaler_ori, clf_ori = Train_ENTRNA_ori()

for row in range(0,len(df)):
    RNA_name = df['Name_0'][row]

    ori_seq = df['Ori_Seq'][row]
    new_seq = Transfer_Ori_Seq_To_Our_Form(ori_seq)
    new_seq = Transfer_T_into_U(new_seq)
    new_seq = ''.join(map(str, new_seq))
    df['Ori_Seq'][row] = new_seq

    if RNA_name != 'NONE':
    #and RNA_name in txt_name_list:
        print 'iteration:', row
        print 'Name:', RNA_name
        #ori_seq = df['Ori_Seq'][row]
        #new_seq = Transfer_Ori_Seq_To_Our_Form(ori_seq)
        #new_seq = Transfer_T_into_U(new_seq)
        #new_seq = ''.join(map(str, new_seq))
    
        Actual_str = df['Actual_str'][row]
        fc = RNA.fold_compound(new_seq)
        (s, mm) = fc.mfe()
        RNAfold_str = s
        RNAfold_str = ModifyingWholeChainByRNAfold(RNAfold_str)
        #df['RNAfold_str'][row] = RNAfold_str
        Our_alg_str = df['Our_alg_str'][row]
    
        Actual_foladability_MFE = Calculate_Foldability(new_seq, Actual_str, scaler_ori, clf_ori)
        Actual_foladability_NFE = Calculate_Foldability_NFE(new_seq, Actual_str, scaler, clf)
        Actual_foladability_abs = abs(Actual_foladability_MFE - Actual_foladability_NFE)
    
        RNAfold_foladability_MFE = Calculate_Foldability(new_seq, RNAfold_str, scaler_ori, clf_ori)
        RNAfold_foladability_NFE = Calculate_Foldability_NFE(new_seq, RNAfold_str, scaler, clf)
        RNAfold_foladability_abs = abs(RNAfold_foladability_MFE - RNAfold_foladability_NFE)
    
        Our_foladability_MFE = Calculate_Foldability(new_seq, Our_alg_str, scaler_ori, clf_ori)
        Our_foladability_NFE = Calculate_Foldability_NFE(new_seq, Our_alg_str, scaler, clf)
        Our_foladability_abs = abs(Our_foladability_MFE - Our_foladability_NFE)
    
        #Changes
        df['Ori_Seq'][row] = new_seq
        df['RNAfold_str'][row] = RNAfold_str
    
        df['Actual_Foldability_MFE'][row] = Actual_foladability_MFE
        df['Actual_Foldability_WOMFE'][row] = Actual_foladability_NFE
        df['Actual_abs_Foldability_difference'][row] = Actual_foladability_abs
    
        df['RNAfold_Foldability_MFE'][row] = RNAfold_foladability_MFE
        df['RNAfold_Foldability_WOMFE'][row] = RNAfold_foladability_NFE
        df['RNAfold_abs_Foldability_difference'][row] = RNAfold_foladability_abs
    
        df['Our_Foldability_MFE'][row] = Our_foladability_MFE
        df['Our_Foldability_WOMFE'][row] = Our_foladability_NFE
        df['Our_abs_Foldability_difference'][row] = Our_foladability_abs
    
df.to_csv(os.path.join(path_of_destination, 'rfam_MFENFE_update_v4_new.csv'))
    
    
    
    