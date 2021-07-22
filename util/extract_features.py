import pandas as pd
import numpy as np
import RNA
import math
from . import rna_toolkit
from .pseudo_decomposition import Peudo_Decom


def extract_features_pseudoknot_free(seq, sec_str):
    feature_dict = {}
    seq_temp = seq
    dp = sec_str
    bp = rna_toolkit.dp_to_bp(dp)
    feature_dict['ent_3'] = []
    feature_dict['gc_percentage'] = []
    feature_dict['ensemble_diversity'] = []
    feature_dict['expected_accuracy'] = []
    #feature_dict['fe_per']  = []

    if (rna_toolkit.is_pseudoknotted(bp) > 0) or (len(bp) != len(seq_temp)):
        feature_dict['ent_3'].append(float('nan'))
        feature_dict['gc_percentage'].append(float('nan'))
        feature_dict['ensemble_diversity'].append(float('nan'))
        feature_dict['expected_accuracy'].append(float('nan'))
        #feature_dict['fe_per'].append(float('nan'))
    else:
        dp_temp = rna_toolkit.bp_to_dp(bp)


        a = RNA.fold_compound(seq_temp)
        a.pf()
        bp_prob = np.array(a.bpp())
        prob_unbp_array = np.ones(len(seq_temp) + 1)
        (s, mm) = a.mfe()
        for ii in range(len(bp) + 1):
            prob_unbp_array[ii] -= np.sum(bp_prob[ii, :])
            prob_unbp_array[ii] -= np.sum(bp_prob[:, ii])

        expected_accuracy = 0
        gamma = 1
        for ii in range(len(seq_temp)):
            if bp[ii] == 0:
                expected_accuracy += prob_unbp_array[ii + 1]
            elif ii + 1 < bp[ii]:
                expected_accuracy += 2 * gamma * bp_prob[ii + 1, bp[ii]]
        expected_accuracy /= len(seq_temp)

        ensemble_diversity = 0
        dim_x, dim_y = bp_prob.shape
        for bp_i in bp_prob[np.triu_indices(dim_x)]:
            ensemble_diversity += 2 * bp_i * (1 - bp_i)
        ensemble_diversity /= len(seq_temp)

        pos_entropy = 0
        for ubp_i in prob_unbp_array:
            if ubp_i > 0 and ubp_i < 1:
                pos_entropy -= (ubp_i * math.log(ubp_i) + (1 - ubp_i) * math.log(1 - ubp_i))
        pos_entropy /= len(seq_temp)

        gc_percentage = 0
        gc_percentage = (seq_temp.count('G') + seq_temp.count('C')) / float(len(seq_temp))

        bp_percentage = 0
        bp_percentage = (dp_temp.count('(') + dp_temp.count(')')) / float(len(seq_temp))

        if rna_toolkit.entropy_max(len(seq_temp), 3) > 0:
            ent_3 = rna_toolkit.entropy(seq_temp, 3) / rna_toolkit.entropy_max(len(seq_temp), 3)
        else:
            ent_3 = float('nan')

        fe_temp = a.eval_structure(dp_temp)
        if mm != 0:
            fe_per = abs(mm - fe_temp) / abs(mm)
        else:
            fe_per = float('nan')

        feature_dict['ent_3'].append(ent_3)
        feature_dict['gc_percentage'].append(gc_percentage)
        feature_dict['ensemble_diversity'].append(ensemble_diversity)
        feature_dict['expected_accuracy'].append(expected_accuracy)
        #feature_dict['fe_per'].append(fe_per)

    df = pd.DataFrame(feature_dict,index=[0])
    return df

def extract_features_pseudoknotted(seq, bp):
    feature_dict = {}
    seq_temp = seq

    feature_dict['gc_percentage'] = float('nan')

    feature_dict['ent_3'] = float('nan')
    feature_dict['ent_4'] = float('nan')
    feature_dict['ent_5'] = float('nan')
    feature_dict['ent_6'] = float('nan')
    feature_dict['ent_7'] = float('nan')
    feature_dict['ent_8'] = float('nan')

    feature_dict['bfe_per']  = float('nan')
    feature_dict['kfe_per']  = float('nan')

    if rna_toolkit.entropy_max(len(seq_temp), 3) > 0:
        ent_3 = rna_toolkit.entropy(seq_temp, 3) / rna_toolkit.entropy_max(len(seq_temp), 3)

    if rna_toolkit.entropy_max(len(seq_temp), 4) > 0:
        ent_4 = rna_toolkit.entropy(seq_temp, 4) / rna_toolkit.entropy_max(len(seq_temp), 4)

    if rna_toolkit.entropy_max(len(seq_temp), 5) > 0:
        ent_5 = rna_toolkit.entropy(seq_temp, 5) / rna_toolkit.entropy_max(len(seq_temp), 5)

    if rna_toolkit.entropy_max(len(seq_temp), 6) > 0:
        ent_6 = rna_toolkit.entropy(seq_temp, 6) / rna_toolkit.entropy_max(len(seq_temp), 6)

    if rna_toolkit.entropy_max(len(seq_temp), 7) > 0:
        ent_7 = rna_toolkit.entropy(seq_temp, 7) / rna_toolkit.entropy_max(len(seq_temp), 7)

    if rna_toolkit.entropy_max(len(seq_temp), 8) > 0:
        ent_8 = rna_toolkit.entropy(seq_temp, 8) / rna_toolkit.entropy_max(len(seq_temp), 8)

    gc_percentage = (seq_temp.count('G') + seq_temp.count('C')) / float(len(seq_temp))
    


    a = RNA.fold_compound(seq_temp)
    a.pf()
    (s, mfe) = a.mfe()

    bfe,kfe = Peudo_Decom(bp,seq,'a')
    if mfe != 0:
        bfe_per = abs(bfe - mfe) / abs(mfe)
        feature_dict['bfe_per']  = bfe_per

    if bfe != 0:
        kfe_per = abs(bfe - kfe) / abs(bfe)
        feature_dict['kfe_per']  = kfe_per


    feature_dict['ent_3'] = ent_3
    feature_dict['ent_4'] = ent_4
    feature_dict['ent_5'] = ent_5
    feature_dict['ent_6'] = ent_6
    feature_dict['ent_7'] = ent_7
    feature_dict['ent_8'] = ent_8

    feature_dict['gc_percentage'] = gc_percentage
    
    df = pd.DataFrame(feature_dict,index=[0])
    return df
