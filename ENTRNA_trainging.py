"""
In this file, we are gonna
(1)Read Seq, and transfer Seq to our form, and Call our alg, and write the our_str into csv (Now hold)
(2)Read Actual str and out them into csv
(3)Read RNAfold str and put them into csv
(4)Calculate distances and put them into csv
(5)Calculate foldability and out them into csv
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
import Main_SeqAttRollOut_RNAfold_3_branches as Our_alg_Main
#*************************************************************

def Read_Name(direct,file_name):
    with open(os.path.join(direct, file_name), "r") as f:
        mylist = f.read().splitlines()
    #print(mylist[0])
    This_name = mylist[0]
    return This_name


def Read_Seq(direct,file_name):
    with open(os.path.join(direct, file_name), "r") as f:
        mylist = f.read().splitlines()
    #print(mylist[1])
    Ori_Seq = mylist[1]
    return Ori_Seq


def Transfer_Ori_Seq_To_Our_Form(Ori_Seq):
    new_seq = [Ori_Seq[i] for i in range(0,len(Ori_Seq))]
    return new_seq





def Read_Actual_str(direct,file_name):
    with open(os.path.join(direct, file_name), "r") as f:
        mylist = f.read().splitlines()
    print(mylist[2][8:])
    Actual_str = mylist[2][8:]
    return Actual_str


def Read_RNAfold_str(direct,file_name):
    with open(os.path.join(direct, file_name), "r") as f:
        mylist = f.read().splitlines()
    print(mylist[3][9:])
    RNAfold_str = mylist[3][9:]
    return RNAfold_str


def Read_Our_alg_str(direct,file_name):
    with open(os.path.join(direct, file_name), "r") as f:
        mylist = f.read().splitlines()
    print(mylist[4][6:])
    RNAfold_str = mylist[4][6:]
    return RNAfold_str


def Calculate_Distance_str(str_actual,str2):
    str_actual = [str_actual[i] for i in range(0, len(str_actual))]
    str2 = [str2[i] for i in range(0, len(str2))]
    #print('str_actual',str_actual)
    #print('str2',str2)
    score = 0.0
    for i in range(len(str_actual)):
        if str_actual[i] != str2[i]:
	    score = score +1
    score = score/float(len(str_actual))
    #print('The distance is:', score)
    return score


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


def Put_something_into_csv(something):
    #python2 can use file instaed of open
    with file("Result_overall_tm_RNA_3_branches.csv", "a") as csvfile: 
        writer = csv.writer(csvfile)
        #writer.writerow([str(something)])
        #writer.writerow(['Name','Ori_Seq', 'Actual_str', 'RNAfold_str', 'Our_alg_str', 'Actual_foladability', 'RNAfold_foladability', 'Our_foldablity', 'RNAfold_distance', 'Our_distance', 'ent_3', 'gc_perentage', 'ensemble_diversity', 'expected_accuracy', 'fe_per'])
        writer.writerow(something)

def Transfer_T_into_U(seq):
    for i in range(len(seq)):
	new_seq = seq
        if new_seq[i] == 'T':
	    new_seq[i] = 'U'
    return new_seq

#**********************************
path_of_ori_data = '/home/menghan/Downloads/ENTRNA-master/Pseudoknot_free/tm_data/'





scaler, clf = Train_ENTRNA()
scaler_ori, clf_ori = Train_ENTRNA_ori()
print('scaler',scaler , 'clf', clf)
print('scaler_ori', scaler_ori, 'clf_ori',clf_ori)







