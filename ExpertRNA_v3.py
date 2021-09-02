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

from util.pseudoknot_free import entrna_main
from util.pseudoknot_free_ori import entrna_main_ori
from util.pseudoknot_free import entrna_train_our_ver
from util.pseudoknot_free_ori import entrna_train_our_ver_ori
from util.pseudoknot_free_ori import entrna_main_return_all_features_ori
import ExpertRNA_main as Our_alg_Main

#*************************************************************

def Read_dbn(directory,file_name):
    with open(os.path.join(directory, file_name), "r") as f:
        lines = f.readlines()
    ext = file_name.split('.')[-1]
    title = file_name.strip(ext)
    seq = lines[1].strip()
    actual = lines[2].strip()
    return (title, seq, actual)

def read_fasta(directory,file_name):
    with open(os.path.join(directory, file_name), "r") as f:
        lines = f.readlines()
    title = lines[0].strip('>').strip()
    seq = lines[1].strip()
    return (title, seq)

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

def testing_header(output_file, expert_nameonly):
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

def testing_output(file_name, Ori_Seq, Actual_str, RNAfold_str, Our_alg_str_list, scaler, clf, scaler_ori, clf_ori, alg_run_time, expert_nameonly):
    #Hamming distance between RNAfold prediction and the actual structure
    RNAfold_distance = Calculate_Distance_str(Actual_str, RNAfold_str)

    # Process the ExpertRNA output and produce output file
    for Our_alg_str in Our_alg_str_list:
        if Our_alg_str != 'NONE':

            #Hamming distance between ExpertRNA prediction and the actual structure
            Our_distance = Calculate_Distance_str(
                Actual_str, Our_alg_str)

            #extract Entrna features of the ExpertRNA prediction
            features_ori = entrna_main_return_all_features_ori(
                Ori_Seq, Our_alg_str, scaler_ori, clf_ori)

            #calculate the free energies of the three structures
            Actual_FE = fc.eval_structure(Actual_str)
            RNAfold_FE = fc.eval_structure(RNAfold_str)
            Expert_FE = fc.eval_structure(Our_alg_str)

            #calculate foldabilities of the three structures
            Actual_foldability_MFE = entrna_main_ori(
                Ori_Seq, Actual_str, scaler_ori, clf_ori)
            Actual_foldability_NFE = entrna_main(
                Ori_Seq, Actual_str, scaler, clf)
            RNAfold_foldability_MFE = entrna_main_ori(
                Ori_Seq, RNAfold_str, scaler_ori, clf_ori)
            RNAfold_foldability_NFE = entrna_main(
                Ori_Seq, RNAfold_str, scaler, clf)
            Expert_foldability_MFE = entrna_main_ori(
                Ori_Seq, Our_alg_str, scaler_ori, clf_ori)
            Expert_foldability_NFE = entrna_main(
                Ori_Seq, Our_alg_str, scaler, clf)

            list_for_the_case = [
                file_name,
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
            list_for_the_case = ['NONE','NONE']
        Put_something_into_csv(list_for_the_case, output_file)
    return

def prod_output(str_name, Ori_Seq, Our_alg_str_list, alg_run_time, output_file, Our_alg_str_list_ori):
    print('Runtime:', alg_run_time)
    print('Writing results to file_name:', output_file)
    with open(output_file, 'a') as f:
        f.write('>'+str_name+'\n')
        f.write(Ori_Seq+'\n')
        for s in Our_alg_str_list_ori:
            f.write(s[0]+'\n')
            f.write(s[1]+'score: '+str(s[2])+'\n')
            

#**********************************
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A RNA structure prediction algorithm which combines multiple heuristics and evaluation functions.')
    parser.add_argument('input_data', type=str, nargs=1, help="The path to a directory containing input sequences.  Each sequence should be a separate fasta file.")
    parser.add_argument('-e', '--expert_set', metavar='expert_set', nargs=2, action='append', help="The Expert used to evaluate partial solutions, and branches to maintain for each expert. Expert should choose from 'ENTRNA_MFE', 'ENTRNA_NFE'. Should in the form of [expert name, branch num] And input should be call the -fd multiple times for each input.")
    parser.add_argument('-f', '--folder_set', metavar='folder_set', nargs=2, action='append', help="The folder to complete partial solutions. Should choose from 'RNAfold' or 'CONTRAfold'. Type should choose from nonspecific or specific. Should in the form of [folder name, type]. And input should be call the -fd multiple times for each input.")
    parser.add_argument('-o', '--output', metavar='output_file', type=str, nargs=1, help='name of file to write results out to in .csv format')
    parser.add_argument('-m', '--min_basepair_distance', metavar='min_bp_dist', type=int, default=4, help='The minimum basepair distance')
    parser.add_argument('-t', '--testing', metavar='testing', action='store_const', const=True, default=False, help='If set the algorithm expects a true structure for each input on the line following the sequence and produces a comparison between the prediction and the real score.')
    args = parser.parse_args()
    
    folder_nameset = args.folder_set
    expert_nameset = args.expert_set
    expert_nameonly = [item[0] for item in expert_nameset]
    
    path_of_ori_data = args.input_data[0]
    min_dbp = args.min_basepair_distance
    
    #parse output file target and number of branches
    if args.output:
        output_file = args.output[0]
    else:
        if args.testing:
            output_file = "ExpertRNA_output_v3_test.csv"
        else:
            output_file = "ExpertRNA_output_v3.dbn"

    #prepare output file
    clear_output(output_file)
    if args.testing:
        testing_header(output_file, expert_nameonly)
    else:
        with open(output_file, 'w+') as f:
            f.write('')

    # Train the ENTRNA model
    scaler, clf = Train_ENTRNA()
    scaler_ori, clf_ori = Train_ENTRNA_ori()



    for file_name in os.listdir(path_of_ori_data):
        try:
            #In this step we tackle data which contains T and substitute T with U
            if args.testing:
                str_name, Ori_Seq, Actual_str = Read_dbn(path_of_ori_data, file_name)
            else:
                str_name, Ori_Seq = read_fasta(path_of_ori_data, file_name)

            new_seq = Transfer_Ori_Seq_To_Our_Form(Ori_Seq)
            new_seq = Transfer_T_into_U(new_seq)
            Ori_Seq = ''.join(map(str, new_seq)) #???
            
            #check
            print('seq:', Ori_Seq)
            print('folder_nameset:', folder_nameset)
            print('expert_nameset:', expert_nameset)
            print('min_dbp:', min_dbp)

            if args.testing:
                #Fold the structure with RNAfold
                fc = RNA.fold_compound(Ori_Seq)
                (s, mm) = fc.mfe()
                RNAfold_str = s
                RNAfold_str = ModifyingWholeChainByRNAfold(RNAfold_str)
                RNAfold_str = ''.join(map(str, RNAfold_str))

            #Fold the structure with ExpertRNA
            alg_start_time = time.time()
            Our_alg_str_list_ori = Our_alg_Main.ExpertRNA(new_seq, folder_nameset, expert_nameset, min_dbp, scaler, clf, scaler_ori, clf_ori)
            #print('Our_alg_str_list_ori:', Our_alg_str_list_ori)
            Our_alg_str_list = [item[0] for item in Our_alg_str_list_ori]
            #print('Our_alg_str_list:', Our_alg_str_list)
            alg_run_time = time.time() - alg_start_time

            if args.testing:
                testing_output(file_name, Ori_Seq, Actual_str, RNAfold_str, Our_alg_str_list, scaler, clf, scaler_ori, clf_ori, alg_run_time, expert_nameonly)
            else:
                prod_output(str_name, Ori_Seq, Our_alg_str_list, alg_run_time, output_file, Our_alg_str_list_ori)
            
        except Exception as e:
            print(e)
            print("Sorry, there's no corresponding action to fortified solution.")
            continue