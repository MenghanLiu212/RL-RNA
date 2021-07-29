"""
In this file, we are gonna
(1)Read Seq, and transfer Seq to our form, and Call our alg, and write the our_str into csv (Now hold)
(2)Read Actual str and out them into csv
(3)Read RNAfold str and put them into csv
(4)Calculate distances and put them into csv
(5)Calculate foldability and out them into csv
"""

import csv
import os
import time
import argparse

import RNA
from numpy import str_

from util.pseudoknot_free import entrna_main
from util.pseudoknot_free_ori import entrna_main_ori
from util.pseudoknot_free import entrna_train_our_ver
from util.pseudoknot_free_ori import entrna_train_our_ver_ori
from util.pseudoknot_free_ori import entrna_main_return_all_features_ori

#*************************************************************

def Read_dbn(directory,file_name):
    with open(os.path.join(directory, file_name), "r") as f:
        lines = f.readlines()
    ext = file_name.split('.')[-1]
    title = file_name.strip(ext)
    seq = lines[1].strip()
    actual = lines[2].strip()
    return (title, seq, actual)

def Transfer_Ori_Seq_To_Our_Form(Ori_Seq):
    return [i for i in Ori_Seq]

def Calculate_Distance_str(str_actual,str2):
    str_actual = [str_actual[i] for i in range(0, len(str_actual))]
    str2 = [str2[i] for i in range(0, len(str2))]
    score = 0.0
    for i in range(len(str_actual)):
        if str_actual[i] != str2[i]:
            score = score +1
    score = score/float(len(str_actual))
    return score

def Train_ENTRNA():
    scaler, clf = entrna_train_our_ver()
    return scaler, clf

def Train_ENTRNA_ori():
    scaler, clf = entrna_train_our_ver_ori()
    return scaler, clf

def clear_output(fname):
    open(fname, 'w').close()

def Put_something_into_csv(something, fname):
    with open(fname, "a") as csvfile:
        writer = csv.writer(csvfile)
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

#**********************************
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A RNA structure prediction algorithm which combines multiple heuristics and evaluation functions.')
    parser.add_argument('input_data', type=str, nargs=1, help="The path to a directory containing input sequences.  Each sequence should be a separate .ct file.")
    parser.add_argument('expert', type=str, nargs=1, help="The Expert used to evaluate partial solutions.  Should be 'Entrna_MFE', 'Entrna_NFE', or 'Entrna_MFENFE'")
    parser.add_argument('-o', '--output', metavar='output_file', type=str, nargs=1, help='name of file to write results out to in .csv format')
    parser.add_argument('-n', '--branches', metavar='n_branches', type=int, default=4, help='The number of Rollout branches to keep at a time')
    args = parser.parse_args()

    # Parse input directory and expert selection
    path_of_ori_data = args.input_data[0]
    if args.expert[0].lower() == 'entrna_mfe':
        import Main_SeqAttRollOut_RNAfold_only_MFE_4_branch_4_min_Agave_cleared as Our_alg_Main
    elif args.expert[0].lower() == 'entrna_nfe':
        import Main_SeqAttRollOut_RNAfold_only_NFE_4_branch_4_min_Agave_cleared as Our_alg_Main
    elif args.expert[0].lower() == 'entrna_mfenfe':
        import Main_SeqAttRollOut_RNAfold_only_MFENFE_2X2_branches_4_min_Agave_cleared as Our_alg_Main
    else:
        print("Unrecognized expert selection.  Please choose either 'Entrna_MFE', 'Entrna_NFE', or 'Entrna_MFENFE'")
        exit()
    
    #parse output file target and number of branches
    if args.output:
        output_file = args.output[0]
    else:
        output_file = "ExpertRNA_output.csv"
    branches = args.branches

    clear_output(output_file)

    start_time = time.time()

    # Train the ENTRNA model
    scaler, clf = Train_ENTRNA()
    scaler_ori, clf_ori = Train_ENTRNA_ori()

    header = [
        'Name',
        'Ori_Seq', 
        'Actual_str', 
        'RNAfold_str', 
        'Our_alg_str', 
        'RNAfold_distance', 
        'Our_distance', 
        'Ent_3', 
        'GC_percentage', 
        'Ensemble_diversity', 
        'Expected_accuracy', 
        'FE_per', 
        'Running_time(sec)',
        'Actual_foldability_MFE',
        'Actual_foldability_NFE',
        'RNAfold_foldability_MFE',
        'RNAfold_foldability_NFE',
        'ExpertRNA_foldability_MFE',
        'ExpertRNA_foldability_NFE',
        'Actual_FE',
        'RNAfold_FE',
        'ExpertRNA_FE'
    ]
    Put_something_into_csv(header, output_file)

    for file_name in os.listdir(path_of_ori_data):
        try:
            #Read and process input data
            str_name, Ori_Seq, Actual_str = Read_dbn(path_of_ori_data, file_name)
            new_seq = Transfer_Ori_Seq_To_Our_Form(Ori_Seq)
            new_seq = Transfer_T_into_U(new_seq)
            Ori_Seq = ''.join(map(str, new_seq)) 

            #read the actual structure
            fc = RNA.fold_compound(Ori_Seq)

            #Fold the structure with RNAfold
            (s, mm) = fc.mfe()
            RNAfold_str = s
            RNAfold_str = ModifyingWholeChainByRNAfold(RNAfold_str)
            RNAfold_str = ''.join(map(str, RNAfold_str))

            #Fold the structure with ExpertRNA
            alg_start_time = time.time()
            Our_alg_str_list = Our_alg_Main.Main_function(new_seq, scaler, clf, scaler_ori, clf_ori, branches)
            alg_run_time = time.time() - alg_start_time

            #Hamming distance between RNAfold prediction and the actual structure
            RNAfold_distance = Calculate_Distance_str(Actual_str, RNAfold_str)

            # Process the ExpertRNA output and produce output file
            for Our_alg_str in Our_alg_str_list:
                if Our_alg_str != 'NONE':

                    #Hamming distance between ExpertRNA prediction and the actual structure
                    Our_distance = Calculate_Distance_str(Actual_str, Our_alg_str)

                    #extract Entrna features of the ExpertRNA prediction
                    features_ori = entrna_main_return_all_features_ori(Ori_Seq, Our_alg_str, scaler_ori, clf_ori)

                    #calculate the free energies of the three structures
                    Actual_FE = fc.eval_structure(Actual_str)
                    RNAfold_FE = fc.eval_structure(RNAfold_str)
                    Expert_FE = fc.eval_structure(Our_alg_str)

                    #calculate foldabilities of the three structures
                    Actual_foldability_MFE = entrna_main_ori(Ori_Seq, Actual_str, scaler_ori, clf_ori)
                    Actual_foldability_NFE = entrna_main(Ori_Seq, Actual_str, scaler, clf)
                    RNAfold_foldability_MFE = entrna_main_ori(Ori_Seq, RNAfold_str, scaler_ori, clf_ori)
                    RNAfold_foldability_NFE = entrna_main(Ori_Seq, RNAfold_str, scaler, clf)
                    Expert_foldability_MFE = entrna_main_ori(Ori_Seq, Our_alg_str, scaler_ori, clf_ori)
                    Expert_foldability_NFE = entrna_main(Ori_Seq, Our_alg_str, scaler, clf)

                    list_for_the_case = [
                        str_name, 
                        Ori_Seq, 
                        Actual_str, 
                        RNAfold_str, 
                        Our_alg_str, 
                        RNAfold_distance, 
                        Our_distance, 
                        features_ori['ent_3'], 
                        features_ori['gc_percentage'], 
                        features_ori['ensemble_diversity'], 
                        features_ori['expected_accuracy'], 
                        features_ori['fe_per'], 
                        alg_run_time,
                        Actual_foldability_MFE,
                        Actual_foldability_NFE,
                        RNAfold_foldability_MFE,
                        RNAfold_foldability_NFE,
                        Expert_foldability_MFE,
                        Expert_foldability_NFE,
                        Actual_FE,
                        RNAfold_FE,
                        Expert_FE
                    ]
                else:
                    list_for_the_case=['NONE']
                Put_something_into_csv(list_for_the_case, output_file)
        except Exception as e:
            print(e)
            print("Sorry, there's no corresponding action to fortified solution.")
            continue
