#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: menghanliu
"""

from copy import deepcopy
from collections import defaultdict
import pandas as pd
import numpy as np
import csv
import json
from pandas import DataFrame
import random
import RNA
import ExpertRNA_toolbox as toolbox
from util.pseudoknot_free import entrna_main
from util.pseudoknot_free_ori import entrna_main_ori
import itertools
import subprocess
import os


from util.pseudoknot_free import entrna_main
from util.pseudoknot_free_ori import entrna_main_ori
from util.pseudoknot_free import entrna_train_our_ver
from util.pseudoknot_free_ori import entrna_train_our_ver_ori
from util.pseudoknot_free_ori import entrna_main_return_all_features_ori

#---------------------------------------------------------------------------------------------------------------------
class Folder():
    def __init__(self, Name, BPInputType):
        """
        For each action, we have four elements of it:
            (1) Name
            (2) BPInputType:can be nonspecific(only open, close or unpair) or specific (with who pair to who
            
        """
        self.Name = Name
        self.BPInputType = BPInputType
        #self.Branch_num = Branch_num
        
    def BaseHeuristics_folder(self, PartialChain, action, OriginalRNAChain, scaler, clf, scaler_ori, clf_ori):
        """
        $$$ NOT IN USE. $$$
        """
        if self.Name.lower() == 'rnafold':
            #RNAfold
            WholeChain = self.BaseHeuristics_RNAfold(self, BPchain, seq, now_position, now_BP, OriginalRNAChain)
        elif self.Name.lower() == 'contrafold':
            #CONTRAfold
            WholeChain = self.BaseHeuristics_CONTRAfold(self, seq, constraint_pair_set, now_position)
        return WholeChain
    
    def BaseHeuristics_RNAfold(self, BPchain, action, seq, now_position, now_BP, OriginalRNAChain, min_dbp):
        """
        In this function, we input our partial chain, which is the full seq and partial BP chain, and we transfer BP chain into hard constraints.
        """
        def ApplyConstraintsFromPartialChain(seq, BPchain):
            """
            In this function, we generate a RNAfold fc model using seq and partial chain BPchain, which contains partial chain info constraints.
            """
            fc = RNA.fold_compound(seq)
            for i in range(len(BPchain)):
                k = i + 1
                if BPchain[i] == 1 or BPchain[i] == 1.5:
                    fc.hc_add_bp_nonspecific(k, 1)
                elif BPchain[i] == 2 or BPchain[i] == 2.5:
                    fc.hc_add_bp_nonspecific(k, -1)
                elif BPchain[i] == 0:
                    fc.hc_add_up(k)
            return fc

        def ModifyingWholeChainByRNAfold(s, now_position):
            for i in range(now_position+1,len(s)-4):
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

        def GetWholeStrChain_withConstrainedRNAfold(seq, BPchain, action):
            """
            In this function, due to iterator problem, we first discuss the length of partial chain.
            And then we get the whole str chain by Constrained RNAfold.
            There are 4 situations we need to consider when adding our constraints to RNAfold in order to disable (...):
            Partial chain ends with ( or (. or (.. or (... ok?
            """
            BPchain_new = deepcopy(BPchain)
            BPchain_new.append(action.BasePair)
            BPchain = BPchain_new
            L = len(BPchain)
            if L == 1:
                if BPchain[-1] == 1:
                    WC = RNAfold_Constrained_Situation(seq, BPchain, L, 1)
                    return WC
                else:
                    fc = ApplyConstraintsFromPartialChain(seq, BPchain)
                    (s, mm) = fc.mfe()
                    return s
            elif L == 2:
                if BPchain[-1] == 1:
                    WC = RNAfold_Constrained_Situation(seq, BPchain, L, 1)
                    return WC
                elif BPchain[-1] == 0 and BPchain[-2] == 1:
                    WC = RNAfold_Constrained_Situation(seq, BPchain, L, 2)
                    return WC
                else:
                    fc = ApplyConstraintsFromPartialChain(seq, BPchain)
                    (s, mm) = fc.mfe()
                    return s
            elif L == 3:
                if BPchain[-1] == 1:
                    WC = RNAfold_Constrained_Situation(seq, BPchain, L, 1)
                    return WC
                elif BPchain[-1] == 0 and BPchain[-2] == 1:
                    WC = RNAfold_Constrained_Situation(seq, BPchain, L, 2)
                    return WC
                elif BPchain[-1] == 0 and BPchain[-2] == 0 and BPchain[-3] == 1:
                    WC = RNAfold_Constrained_Situation(seq, BPchain, L, 3)
                    return WC
                else:
                    fc = ApplyConstraintsFromPartialChain(seq, BPchain)
                    (s, mm) = fc.mfe()
                    return s
            elif L >= 4:
                if BPchain[-1] == 1:
                    WC = RNAfold_Constrained_Situation(seq, BPchain, L, 1)
                    return WC
                elif BPchain[-1] == 0 and BPchain[-2] == 1:
                    WC = RNAfold_Constrained_Situation(seq, BPchain, L, 2)
                    return WC
                elif BPchain[-1] == 0 and BPchain[-2] == 0 and BPchain[-3] == 1:
                    WC = RNAfold_Constrained_Situation(seq, BPchain, L, 3)
                    return WC
                elif BPchain[-1] == 0 and BPchain[-2] == 0 and BPchain[-3] == 0 and BPchain[-4] == 1:
                    WC = RNAfold_Constrained_Situation(seq, BPchain, L, 4)
                    return WC
                else:
                    fc = ApplyConstraintsFromPartialChain(seq, BPchain)
                    (s, mm) = fc.mfe()
                    return s    
            else:
                raise Exception("Sorry, should never reach here, the length of BPchain is wrong")

        def RNAfold_Constrained_Situation(seq, BPchain, L, situation_num):
            """
            In this function, we generate all possible combination of base-pair in situation num, and select the whole chain str with min MFE.
            In situation 1, there is a ( at the end of our partial chain, so there are four positions should not be ).
            So there are 2^4 possible combinations of . and ( at these four positions.
            So similar for other situation_num.
            We generate all of them, and argmin_WC (FreeEnergy).
            """
            if situation_num == 1:
                n = 4
            elif situation_num == 2:
                n = 3
            elif situation_num == 3:
                n = 2
            elif situation_num == 4:
                n = 1
            else:
                raise Exception("Sorry, should never reach here!")
            WC_set = []
            fc_set = ApplyConstraintsForAllPossibleCombo(seq, n, L, BPchain)
            for item in fc_set:
                (s,m) = item.mfe()
                FreeEnergy = item.eval_structure(s)
                WC_set.append([s, FreeEnergy])
            removal_list = []
            for item in WC_set:
                dec = DecideIfAllDot(item[0])
                if dec == True:
                    removal_list.append(item)
            if removal_list != []:
                for item in removal_list:
                    WC_set.remove(item)
            if WC_set == []:
                print('Empty WC_set! No feasible result from RNAfold for this partial chain and action!')
                best_WC = 'No feasible result from RNAfold'
            else:
                best_WC = min(WC_set, key= lambda x: x[1])[0]
            return best_WC

        def DecideIfAllDot(s):
            temp=0
            for i in range(0,len(s)):
                if s[i] !='.':
                    temp = temp+1
                else:
                    pass
            if temp ==0:
                decision = True
            else:
                decision = False
            return decision

        def ApplyConstraintsForAllPossibleCombo(Ori_seq, num_of_postions, L, BPchain):
            """
            In this function, we (1) generate code represents the arrangement of constraints, (2) translate the codes into constraints.
            """
            fc_set = []
            code_list = list(itertools.product([0,1], repeat=num_of_postions))
            for item in code_list:
                seq_temp = deepcopy(Ori_seq)
                fc_temp = ApplyConstraintsFromPartialChain(seq_temp, BPchain)
                for i in range(len(item)):
                    if item[i] == 0:
                        fc_temp.hc_add_up(L + 1 + i)
                    elif item[i] == 1:
                        fc_temp.hc_add_bp_nonspecific(L + 1 + i, 1)
                fc_set.append(fc_temp)
            return fc_set

        #function begin
        print('--*--BaseHeuristics_RNAfold--*--')
        if min_dbp == 4:
            #start restricting for min_dbp =4
            s = GetWholeStrChain_withConstrainedRNAfold(seq, BPchain, action)
            if s != 'No feasible result from RNAfold':
                new_s = [s[i] for i in range(0, len(s))]
                new_s = ModifyingWholeChainByRNAfold(new_s, now_position)
                s = ''.join(map(str, new_s))
                #print('seq and str:', seq, s)
            else:
                s='No feasible result from RNAfold'
        else:
            #if min_dbp<=3
            fc = ApplyConstraintsFromPartialChain(seq, BPchain)
            (s, mm) = fc.mfe()
            #print("%s\n%s (MM: %d)\n" %  (seq, s, -mm))
        return s


    def BaseHeuristics_CONTRAfold(self, seq, constraint_pair_set, now_position):
        # In this function, this is the direct CONTRAfold application, input action with who pair to who, and solution
        #This is the function of CONTRAfold,  we need to turn pc and constraints into BPSEQ file and save it to a directory, then extract
        #print('--*--BaseHeuristics_CONTRAfold--*--')
        #print('constraint_pair_set:', constraint_pair_set)
        #print('now_position:', now_position)
        #print 'seq:', seq
        file_path = "/home/mliu126/RNA/RL-RNA/Alg_running/3_branches/ENTRNA-master_new_MH/ENTRNA-master/CONTRAfold_need_temp_saving/temp_saving_CONTRAfold/"
        file_name = "generated_RNA.BPSEQ"

        with open(file_path + file_name, "a") as txtfile:
        #with file(file_path + file_name, "a") as txtfile: #python2 ver
            for row in range(0, len(seq)):
                txtfile.write(str(row+1))
                txtfile.write(' ')
                txtfile.write(seq[row])
                txtfile.write(' ')
                temp = 0
                for item in constraint_pair_set:
                    if row == item[0]:
                        res = item[1]+1
                        temp += 1
                    elif row == item[1]:
                        res = item[0]+1
                        temp += 1
                    else:
                        pass
                if temp == 0:
                    if row <= now_position:
                        res = 0
                    else:
                        res = -1
                txtfile.write(str(res))
                txtfile.write("\n")

        our_file = file_path + file_name
        our_command = "./contrafold predict " + our_file + " --constraints"
        #output = subprocess.Popen(our_command, stdout=subprocess.PIPE, shell=True, cwd="/home/mliu126/contrafold-se/src").communicate() #for agave
        output = subprocess.Popen(our_command, stdout=subprocess.PIPE, shell=True, cwd="/home/menghan/Downloads/contrafold/src").communicate()  #for pc
        #print('output:', output)
        output_split_list = output[0].splitlines()
        #WC_dp = output_split_list[2]
        WC_dp = output_split_list[2]
        #position 2 for pc and position 4 for Agave
        #print('structure:', WC)

        #delete this file
        os.remove(our_file)
        return WC_dp

        
class Expert():
    def __init__(self, Name, Branch_num):
        """
        For each action, we have four elements of it:
            (1) Name
            (2) Branch_num: num of branches assigned to this expert
        """
        self.Name = Name
        #self.Type = Type
        self.Branch_num = Branch_num

    def GetRewards(self, dp_str, OriginalRNAChain, scaler, clf, scaler_ori, clf_ori):
        seq = OriginalRNAChain
        seq_str = ''.join(map(str, seq))
        if self.Name.lower() == 'entrna_mfe':
            foldability = entrna_main_ori(seq_str, dp_str, scaler_ori, clf_ori)
            rewards = foldability
        elif self.Name.lower() == 'entrna_nfe':
            foldability = entrna_main(seq_str, dp_str, scaler, clf)
            rewards = foldability
        else:
            #more experts to be included
            pass
        return rewards


#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
def UpdateSolutionSet(itera, SolutionSet, RNAOriginalChain, Folder_set, Expert_set, min_dbp, scaler, clf, scaler_ori, clf_ori):
    print('--*--UpdateSolutionSet--*--')
    GlobalPossibleActionSet = []
    for pc in SolutionSet:
        LocalPossibleActionSet = GenerateLocalPossibleActionSet(pc, RNAOriginalChain, min_dbp)
        #update rewards for each local action set
        LocalPossibleActionSet = UpdateRewardsForAction(LocalPossibleActionSet, RNAOriginalChain, Folder_set, Expert_set, min_dbp, scaler, clf, scaler_ori, clf_ori)
        pc = UpdateFortifiedSolutionForEachSolution(pc, LocalPossibleActionSet, Folder_set, Expert_set)
        print('*LocalPossibleActionSet:*')
        for lpa in LocalPossibleActionSet:
            print('Action:', 'Nucleotide:', lpa.Nucleotide, 'BasePair:', lpa.BasePair, 'PositionOfNucleotide:',lpa.PositionOfNucleotide, 'Reward', lpa.rewards, 'PartialChain_parent_BP_Chain', lpa.PartialChain_parent.BasePairChain)
            GlobalPossibleActionSet.append(lpa)
    print('*GlobalPossibleActionSet:*')
    for gpa in GlobalPossibleActionSet:
        print('Action:', 'Nucleotide:', gpa.Nucleotide, 'BasePair:', gpa.BasePair, 'PositionOfNucleotide:',
              gpa.PositionOfNucleotide, 'Reward', gpa.rewards, 'PartialChain_parent_BP_Chain', gpa.PartialChain_parent.BasePairChain)

    TopActionSet = GenereateTopActionSet(GlobalPossibleActionSet, SolutionSet, RNAOriginalChain, Folder_set, Expert_set, scaler, clf, scaler_ori, clf_ori)
    print('*TopActionSet*')
    for ta_full in TopActionSet:
        ta=ta_full[0]
        print('Action:', 'Nucleotide:', ta.Nucleotide, 'BasePair:', ta.BasePair,
        'rewards', ta.rewards, 'PartialChain_parent_BP_Chain', ta.PartialChain_parent.BasePairChain)

    new_SolutionSet = []
    for item in TopActionSet:
        action = item[0]
        this_ept = item[1]
        this_ept_score = item[2]
        New_PartialChain = toolbox.TakeAction(action.PartialChain_parent, action)
        new_SolutionSet.append([New_PartialChain, this_ept, this_ept_score])
    return new_SolutionSet


def UpdateFortifiedSolutionForEachSolution(Solution, LocalPossibleActionSet, Folder_set, Expert_set):
    """
    Update the Fortified solution matrix for pc according to action and WC.
    """
    if Solution.ori_WC == 'NNN':
        ori_WC_dict = {}
        action_init = random.choice(LocalPossibleActionSet)
        #print('this test one:', action_init.ori_WC)
        for ept in Expert_set:
            for fd in Folder_set:
                #ori_WC_dict.setdefault(ept.Name, {})[fd.Name] = [action_init.ori_WC[ept.Name][fd.Name][0], action_init.ori_WC[ept.Name][fd.Name][1]]
                #print('this test one:', action_init.ori_WC)
                ori_WC_dict.setdefault(ept.Name, {})[fd.Name] = action_init.rewards[ept.Name][fd.Name]
        #print('ori_WC_dict:', ori_WC_dict)
        Solution.ori_WC = ori_WC_dict
    else:
        pass
    for action in LocalPossibleActionSet:
        for ept in Expert_set:
            for fd in Folder_set:
                #print('Solution.ori_WC:', Solution.ori_WC)
                old = Solution.ori_WC[ept.Name][fd.Name]
                this = action.rewards[ept.Name][fd.Name]
                if this[1] >= old[1]:
                    Solution.ori_WC[ept.Name][fd.Name] = this
        
    return Solution


def GenerateLocalPossibleActionSet(Solution, RNAOriginalChain, min_dbp):
    """
    For each PartialChain in the PartialChainSet, we generate its feasible actions, and also keep the trace(memo) of the PartialChain itself.
    So we get a LocalPossibleActionSet, which is an attribute of the ParitalChain.
    Also, we append all the actions in LocalPossibleActionSet to GlobalPossibleActionSet.
    """
    print('--*--GenerateLocalPossibleActionSet--*--')
    LocalPossibleActionSet = toolbox.GetPossibleActions(Solution, RNAOriginalChain, min_dbp)
    return LocalPossibleActionSet


def UpdateRewardsForAction(PossibleActionSet, OriginalRNAChain, Folder_set, Expert_set, min_dbp, scaler, clf, scaler_ori, clf_ori):
    print('--*--UpdateRewardsForAction--*--')
    PossibleActionSet_new=[]
    for action in PossibleActionSet:
        Solution = action.PartialChain_parent
        action = RollOut(Solution, action, OriginalRNAChain, Folder_set, Expert_set, min_dbp, scaler, clf, scaler_ori, clf_ori)
        if action.ori_WC == 'RNAfold infeasible':
            print('***** RNAfold infeasible, cut! *****')
            PossibleActionSet.remove(action)
        elif action.ori_WC == 'CONTRAfold infeasible':
            print('***** CONTRAfold infeasible, cut! *****')
            PossibleActionSet.remove(action)
        else:
            #calculate rewards for each expert for each WC from each folder, for each action
            #the reward is in the form of {'ExpertName': {'FolderName':[WC, reward]}}
            rwd_set = {}
            for ept in Expert_set:
                for fd_name, WC in action.ori_WC.items():
                    reward = ept.GetRewards(WC, OriginalRNAChain, scaler, clf, scaler_ori, clf_ori)
                    rwd_set.setdefault(ept.Name, {})[fd_name] = [WC,reward]
                action.rewards = rwd_set
        #print('updated action rewards:', action.rewards)
        PossibleActionSet_new.append(action)
    return PossibleActionSet_new



def RollOut(Solution, action, OriginalRNAChain, Folder_set, Expert_set, min_dbp, scaler, clf, scaler_ori, clf_ori):
    """
    In here, for each action, we do rollout by each folder, and finally select the best wc as a rollout result.
    Also for the CONTRAfold part, we need to decide if this action is feasible for CONTRAfold
    """
    def CONTRAfold_Permute_Selection(PC, seq, WC_before, now_position, now_BP, OriginalRNAChain, scaler, clf, scaler_ori, clf_ori, PartialChainSim, min_dbp):
        # In this function, we permuatate all possible pairing for this action pos, and select the best pairing and WC(CONTRAfold) by ENTRNA
        #print('--*--CONTRAfold_Permute_Selection:--*--')
        #PartialChainCopy = deepcopy(SolutionSim)
        WC_before_bp = toolbox.TransferOursIntoDP(WC_before)
        PCPairingInfo = toolbox.TurnWholeChainToPCPairingInfo(WC_before_bp, PC)
        print('PCPairingInfo before rollout:', PCPairingInfo)
        if now_BP == 0:
            chosen_pair = None
            chosen_WC = fd.BaseHeuristics_CONTRAfold(seq, PCPairingInfo, now_position)
        elif now_BP == 2:
            open_pair = len(PC)-1-PC[::-1].index(1)
            #print('open pair:', open_pair)
            #delete the pair with last open
            remove_list=[]
            for item in PCPairingInfo:
                print(item)
                if item[0]==open_pair or item[1]==now_position:
                    print('remove item')
                    remove_list.append(item)
            for item in remove_list:
                PCPairingInfo.remove(item)
            chosen_pair = [open_pair, now_position]
            PCPairingInfo.append(chosen_pair)
            chosen_WC = fd.BaseHeuristics_CONTRAfold(seq, PCPairingInfo, now_position)
        elif now_BP == 1:
            if not PCPairingInfo:
                permutate_range = [now_position+5, len(seq)]
            else:
                nearest_pair = max(PCPairingInfo, key = lambda i: i[0])
                permutate_range = [now_position+1+min_dbp, nearest_pair[1]]
            choice_set=[] #should be composed of [this_pair, this_WC, rewards]
            print('permutate_range:', permutate_range[0], permutate_range[1], 'length of OriginalRNAChain:', len(OriginalRNAChain))
            feasible_pos_set = []
            for p in range(permutate_range[0], permutate_range[1]):
                check = toolbox.CheckIfNucleotidesCanBePaired(now_position, p, OriginalRNAChain)
                if check == 1:
                    feasible_pos_set.append(p)
            print('feasible_pos_set:', feasible_pos_set)
            for pos in feasible_pos_set:
                this_pair = [now_position, pos]
                PCPairingInfo_new = deepcopy(PCPairingInfo)
                PCPairingInfo_new.append(this_pair)
                this_WC = fd.BaseHeuristics_CONTRAfold(seq, PCPairingInfo_new, now_position)
                #rewards_foldability_withoutMFE = GetRewards(this_WC, OriginalRNAChain, scaler, clf)
                rewards_foldability_withMFE = GetRewards_ori(this_WC, OriginalRNAChain, scaler_ori, clf_ori)
                choice_set.append([this_pair, this_WC, rewards_foldability_withMFE])
            chosen_pair = max(choice_set, key = lambda i: i[2])[0]
            chosen_WC = max(choice_set, key = lambda i: i[2])[1]
            return chosen_pair, chosen_WC
    
    def GetRewards_ori(dp_str, OriginalRNAChain, scaler_ori, clf_ori):
        """
        In this function we plug in ENTRNA as rewards.
        $$$TEMPORARY$$$
        """
        seq = OriginalRNAChain
        #dp = SeqAttRollOut_RNAfold.TransferIntoDP(WholeChain)
        seq_str = ''.join(map(str, seq))
        #dp_str = ''.join(map(str, dp))
        foldability = entrna_main_ori(seq_str, dp_str, scaler_ori, clf_ori)
        rewards = foldability
        #print('Now the rewards is',rewards)
        return rewards

    #Rollout function begin
    print('--*--Rollout:--*--')
    SolutionSim = deepcopy(Solution)
    updated_action = deepcopy(action)
    now_position = action.PositionOfNucleotide
    now_BP = action.BasePair
    # seq = PartialChainSim.ChainItself
    BPchain = SolutionSim.BasePairChain
    WC_before = SolutionSim.ori_WC
    seq = ''.join(map(str, OriginalRNAChain))
    
    for fd in Folder_set:
        if fd.Name.lower() == 'rnafold':
            s = fd.BaseHeuristics_RNAfold(BPchain, action, seq, now_position, now_BP, OriginalRNAChain, min_dbp)
            #check if the chain s is an unfolded chain, if yes, we use fortified_WC from the partial chain
            if s == 'No feasible result from RNAfold':
                WC_dp = 'RNAfold infeasible'
                print('******RNAfold infeasible*****')
                #Solution.ori_WC = 'RNAfold infeasible'
            else:
                WC_dp = s
        elif fd.Name.lower() == 'contrafold':
            #first decide if this action is feasible for CONTRAfold
            if CONTRAfold_feasibility_decision(Solution, action, OriginalRNAChain) == True:
                # Now generate the set of PossibleCombo according to the action and partial chain 
                # And Select the PossibleCombo with highest folder rewards, in here, for RNAfold, is the lowest free energy.
                chosen_pair, chosen_WC = CONTRAfold_Permute_Selection(BPchain, seq, WC_before, now_position, now_BP, OriginalRNAChain, scaler, clf, scaler_ori, clf_ori, SolutionSim)
                s= chosen_WC
                #print('***JUST CHECKING***')
                #print(seq, s)
                #check if the chain s is an unfolded chain, if yes, we use fortified_WC from the partial chain
                if toolbox.DecideIfAllDot(s) == True:
                    WC_dp = 'CONTRAfold infeasible'
                    print('******CONTRAfold infeasible*****')
                else:
                    WC_dp = s
        else:
            #Can add more folders here in the future.
            pass
        updated_action.ori_WC[fd.Name] = WC_dp
        
    return updated_action


def CONTRAfold_feasibility_decision(Solution, action, OriginalRNAChain):
    #decide if this action is feasible for CONTRAfold
    #This part is only for CONTRAfold
    action_BPchoice = action.BasePair
    CONTRAfold_feasibility = True
    WC_before = Solution.ori_WC
    p = len(Solution.ChainItself)
    if action_BPchoice == 1:
        #print 'flag bpcl 1'
        now_position = p
        seq = ''.join(map(str, OriginalRNAChain))
        #WC_before = PartialChain.ori_WC
        WC_before_bp = toolbox.TransferOursIntoDP(WC_before)
        PCPairingInfo = toolbox.TurnWholeChainToPCPairingInfo(WC_before_bp, Solution.BasePairChain)
        #print 'PCPairingInfo:', PCPairingInfo
        if not PCPairingInfo:
            permutate_range = [now_position+5, len(seq)]
        else:
            nearest_pair = max(PCPairingInfo, key = lambda i: i[0])
            permutate_range = [now_position+5, nearest_pair[1]]
        choice_set=[] #should be composed of [this_pair, this_WC, rewards]
        #print('permutate_range:', permutate_range[0], permutate_range[1], 'length of OriginalRNAChain:', len(OriginalRNAChain))
        feasible_pos_set = []
        for p in range(permutate_range[0], permutate_range[1]):
            check = toolbox.CheckIfNucleotidesCanBePaired(now_position, p, OriginalRNAChain)
            if check == 1:
                feasible_pos_set.append(p)
        #print 'feasible_pos_set:', feasible_pos_set
        if feasible_pos_set == []:
            CONTRAfold_feasibility = False
    elif action_BPchoice == 2:
        #conditions for close considering known bp info
        PC_clip = WC_before[0:p]
        if 1 not in PC_clip:
            CONTRAfold_feasibility = False
        pass
    else:
        pass
    return CONTRAfold_feasibility


def GenereateTopActionSet(GlobalPossibleActionSet, SolutionSet, OriginalRNAChain, Folder_set, Expert_set, scaler, clf, scaler_ori, clf_ori):
    """
    We generate the TopActionSet, which is all the top actions with respect to rank in each critics.
    We don't care which action is chosen by which critic.

    We create a temp_WC_set, which is like [[WC,[reward1,reward2,reward3,...]],[WC,[reward1,reward2,reward3,...]]]
    First we generate fortified action from the fortified WC info from PC and add them to set.
    Then we compare and select top n.
    Each expert select its own assigned selection num

    """
    br_num_total = 0
    for ept in Expert_set:
        br_num_total += int(ept.Branch_num)
    #print('Total branch num:', br_num_total)
    
    temp_GlobalPossibleActionSet = []
    for item in GlobalPossibleActionSet:
        temp_GlobalPossibleActionSet.append(item)

    #Consider fortified
    #Make fortified chains into action format and put them into the same pond as other actions
    for PC in SolutionSet:
        if len(PC.BasePairChain) >0 and len(PC.BasePairChain) < len(OriginalRNAChain):    #to avoid the initial chain being NNN
            now_pos = len(PC.BasePairChain)
            try:
                BP_choice = PC.ori_WC[Expert_set[0].Name][Folder_set[0].Name][0][now_pos]
            except:
                print("Dead branch, moving on...")
            if BP_choice == '(':
                BP_choice_bp =1
            elif BP_choice == '.':
                BP_choice_bp =0
            elif BP_choice == ')':
                BP_choice_bp =2
            action_fortified = toolbox.Action(Nucleotide=OriginalRNAChain[now_pos], BasePair=BP_choice_bp,
                                  PositionOfNucleotide=now_pos, PartialChain_parent = PC)
            action_fortified.rewards = PC.ori_WC
            for fd in Folder_set:
                action_fortified.ori_WC[fd.Name] = PC.ori_WC[Expert_set[0].Name][fd.Name][0]
            temp_GlobalPossibleActionSet.append(action_fortified)

    #delete repeated elements
    res = []
    for action in temp_GlobalPossibleActionSet:
        ct = 0
        for item in res:
            if action.Nucleotide == item.Nucleotide and action.BasePair==item.BasePair and action.PositionOfNucleotide==item.PositionOfNucleotide and action.rewards==item.rewards and action.PartialChain_parent.ori_WC==item.PartialChain_parent.ori_WC:
                ct+=1
        if ct == 0:
            res.append(action)
    temp_GlobalPossibleActionSet = res
    #print('temp_GlobalPossibleActionSet after modification:', temp_GlobalPossibleActionSet)

    TopActionSet_temp = []
    if len(temp_GlobalPossibleActionSet) < br_num_total:
        #TopActionSet_temp = temp_GlobalPossibleActionSet
        for item in temp_GlobalPossibleActionSet:
            TopActionSet_temp.append([item, 'not important', 'not important'])
    else:
        UpdatedExpert_set = random.sample(Expert_set, len(Expert_set))
        for ept in UpdatedExpert_set:
            n = int(ept.Branch_num)
            for itera in range(0,n):
                best_res = SelectBestAction(temp_GlobalPossibleActionSet, ept, Folder_set)
                BestForRewards_top_temp = best_res[0]
                best_score = best_res[1]
                temp_GlobalPossibleActionSet.remove(BestForRewards_top_temp)
                TopActionSet_temp.append([BestForRewards_top_temp, ept.Name, best_score])
    return TopActionSet_temp


def SelectBestAction(PossibleActionSet, ept, Folder_set):
    """
    In this fuction we choose the action with best reward.
    """
    temp = PossibleActionSet[0]
    for item in PossibleActionSet:
        reward = max([item.rewards[ept.Name][fd.Name][1] for fd in Folder_set])
        reward_temp = max([temp.rewards[ept.Name][fd.Name][1] for fd in Folder_set])
        if reward > reward_temp:
            temp = item
            reward_temp = reward
    return [temp, reward_temp]


def WholeChainElementToAction(element):
    if element == 1.5:
        res = 1
    elif element == 2.5:
        res = 2
    else:
        res = element
    return res



# -------------------------------------------------------------------------------------
def ExpertRNA(OriginalRNAChain, folder_nameset, expert_nameset, min_dbp, scaler, clf, scaler_ori, clf_ori):
    Folder_set = []
    Expert_set = []
    for item in folder_nameset:
        Folder_set.append(Folder(Name=item[0], BPInputType=item[1]))
    for item in expert_nameset:
        Expert_set.append(Expert(Name=item[0], Branch_num=item[1]))
    
    RNAOriginalChain = OriginalRNAChain
    iteration_length = len(RNAOriginalChain)
    SolutionSet = []

    InitialSolution = toolbox.Solution(ChainItself=[], BasePairChain=[])
    print('***InitialPartialChainSolution***', InitialSolution.ChainItself, InitialSolution.BasePairChain)
    SolutionSet.append(InitialSolution)

    for itera in range(iteration_length):
        print('**************************Begin Iteration*************************************')
        print('**Iteration**', itera)
        print('**OriginalRNAChain**', OriginalRNAChain)
        SolutionSet_with_ept_into = UpdateSolutionSet(itera, SolutionSet, RNAOriginalChain, Folder_set, Expert_set, min_dbp, scaler, clf, scaler_ori, clf_ori)
        SolutionSet = [item[0] for item in SolutionSet_with_ept_into]
        print('***ParitalChainSet:***')
        for item in SolutionSet:
            print('PartialChain:', 'BP_Chain:', item.BasePairChain, 'DP_Chain', toolbox.BP_to_DP(item.BasePairChain), 'Whole chain:', item.ori_WC)
    print('**************************END OF ITER*****************************************')

    #print('**************************Statistics*****************************************')
    #print('the length of RNAchain:', len(RNAOriginalChain))
    Our_dp_str_list = []
    #print('Foldability:')
    for item in SolutionSet_with_ept_into:
        sol = item[0]
        #print('Foldability for this chain:', GetRewards_ori(toolbox.BP_to_DP(item.ori_WC), OriginalRNAChain, scaler_ori, clf_ori))
        Our_BP = sol.BasePairChain
        Our_BPCon = ''.join(map(str, Our_BP))
        Our_dp = toolbox.BP_to_DP(Our_BP)
        Our_dp_str = ''.join(map(str, Our_dp))
        Our_ept = item[1]
        Our_score = item[2]
        Our_dp_str_list.append([Our_dp_str, Our_ept, Our_score])
    #print('Our_dp_str_list:', Our_dp_str_list)
    print('*******************************************************************')
    print('*******************************************************************')
    return Our_dp_str_list


def Train_ENTRNA():
    scaler, clf = entrna_train_our_ver()
    return scaler, clf

def Train_ENTRNA_ori():
    scaler, clf = entrna_train_our_ver_ori()
    return scaler, clf


#test running of the alg
"""
OriginalRNAChain_old= 'UAGGUUUGGCGGUCAUAGCGAUGGGGUUACACCUGGUCUCGUUUCGAUCCCAGAAGUUAAGUCUUUUCGCGUUUUGUUUGUGUACUAUGGGUUCCGGUCUAUGGGAAUUUCAUUUAGCUGCCAGCUUUUU'
OriginalRNAChain1 = [i for i in OriginalRNAChain_old]
folder_nameset= [['RNAfold', 'nonspecific']]
expert_nameset= [['ENTRNA_NFE', 4]]
min_dbp=4
# Train the ENTRNA model
scaler, clf = Train_ENTRNA()
scaler_ori, clf_ori = Train_ENTRNA_ori()

our_res = ExpertRNA(OriginalRNAChain1, folder_nameset, expert_nameset, min_dbp, scaler, clf, scaler_ori, clf_ori)
print('Our result:', our_res)
"""