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
import csv
import os
import Main_SeqAttRollOut_RNAfold_only_MFE_4_branch_4_min_Agave_cleared as Our_alg_Main
import time
import argparse
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
    return [i for i in Ori_Seq]

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


def Calculate_Foldability_NFE(seq, dp_str, scaler, clf):
    print('*Calculate Foldability NFE')
    print('seq:',seq)
    print('dp_str:', dp_str)
    foldability_WOMFE = entrna_main(seq, dp_str, scaler, clf)
    #foldability = entrna_main(seq1,dp_str,scaler, clf)
    #print('foldability',foldability)
    return foldability_WOMFE

def clear_output(fname):
    open(fname, 'w').close()

def Put_something_into_txt(title, something, fname):
    with open(fname, "a") as txtfile:
        txtfile.write(title)
        txtfile.write("\n") 
        txtfile.write(something)
        txtfile.write("\n")

def Put_something_into_csv(something, fname):
    with open(fname, "a") as csvfile:
        writer = csv.writer(csvfile)
        #writer.writerow([str(something)])
        #writer.writerow(['Name','Ori_Seq', 'Actual_str', 'RNAfold_str', 'Our_alg_str', 'Actual_foladability', 'RNAfold_foladability', 'Our_foldablity', 'RNAfold_distance', 'Our_distance', 'ent_3', 'gc_perentage', 'ensemble_diversity', 'expected_accuracy', 'fe_per'])
        writer.writerow(something)

def Transfer_T_into_U(seq):
    return ['U' if i == 'T' else i for i in seq]

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
    return s


def FreeEnergyCalculation(seq, str):
    fc = RNA.fold_compound(seq)
    FE = fc.eval_structure(str)
    return FE


#**********************************
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A RNA structure prediction algorithm which combines multiple heuristics and evaluation functions.')
    parser.add_argument('input_data', type=str, nargs=1, help="The path to a directory containing input sequences.  Each sequence should be a separate .ct file.")
    parser.add_argument('-o', '--output', metavar='output_file', type=str, nargs=1, help='name of file to write results out to in .csv format')
    parser.add_argument('-n', '--branches', metavar='n_branches', type=int, default=4, help='The number of Rollout branches to keep at a time')
    args = parser.parse_args()

    path_of_ori_data = args.input_data[0]
    if args.output:
        output_file = args.output[0]
    else:
        output_file = "ExpertRNA_output.csv"
    branches = args.branches
    #path_of_ori_data = "./Testing_data/5s_RNA_data"

    clear_output(output_file)

    start_time = time.time()

    scaler, clf = Train_ENTRNA()
    scaler_ori, clf_ori = Train_ENTRNA_ori()
    #print('scaler, clf',scaler, clf)
    #print('scaler_ori, clf_ori',scaler_ori, clf_ori)

    title = ['Name','Ori_Seq', 'Actual_str', 'RNAfold_str', 'Our_alg_str', 'RNAfold_distance', 'Our_distance', 'ent_3', 'gc_perentage', 'ensemble_diversity', 'expected_accuracy', 'fe_per', 'Running_time(sec)']
    Put_something_into_csv(title, output_file)

    for file_name in os.listdir(path_of_ori_data):
        
            #scaler, clf = Train_ENTRNA()
            #The_name = Read_Name(path_of_ori_data,file_name)
            The_name = file_name
            #In this step we taclke data which contains T and substitute T with U
            Ori_Seq = Read_Seq(path_of_ori_data,file_name)
            new_seq = Transfer_Ori_Seq_To_Our_Form(Ori_Seq)
            new_seq = Transfer_T_into_U(new_seq)
            Ori_Seq = ''.join(map(str, new_seq)) #???

            Actual_str = Read_Actual_str(path_of_ori_data, file_name)
            fc = RNA.fold_compound(Ori_Seq)
            (s, mm) = fc.mfe()
            RNAfold_str = s
            RNAfold_str = ModifyingWholeChainByRNAfold(RNAfold_str)
            RNAfold_str = ''.join(map(str, RNAfold_str))
            #RNAfold_str = Read_RNAfold_str(path_of_ori_data,file_name)
            #RNAfold_str = ModifyingWholeChainByRNAfold(RNAfold_str)

            alg_start_time = time.time()
            Our_alg_str_list = Our_alg_Main.Main_function(new_seq, scaler, clf, scaler_ori, clf_ori, branches)
            alg_run_time = time.time() - alg_start_time

            #Actual_foladability_MFE = Calculate_Foldability(Ori_Seq,Actual_str, scaler_ori, clf_ori)

            #RNAfold_foladability_MFE = Calculate_Foldability(Ori_Seq, RNAfold_str, scaler_ori, clf_ori)
            RNAfold_distance = Calculate_Distance_str(Actual_str, RNAfold_str)

            for Our_alg_str in Our_alg_str_list:
                if Our_alg_str != 'NONE':
                    #Our_foladability_MFE = Calculate_Foldability(Ori_Seq, Our_alg_str, scaler_ori, clf_ori)

                    Our_distance = Calculate_Distance_str(Actual_str, Our_alg_str)

                    features_ori = entrna_main_return_all_features_ori(Ori_Seq, Our_alg_str, scaler_ori, clf_ori)
                    print(features_ori)

                    list_for_the_case=[The_name, Ori_Seq, Actual_str, RNAfold_str, Our_alg_str, RNAfold_distance, Our_distance, features_ori[0], features_ori[1], features_ori[2], features_ori[3], features_ori[4], alg_run_time]
                else:
                    list_for_the_case='NONE'
                Put_something_into_csv(list_for_the_case, output_file)
        #except Exception as e:
        #    print("Sorry, there's no corresponding action to fortified solution.")
        #    continue







