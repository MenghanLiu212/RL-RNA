
"""
3 branches roll-out alg

We keep all feasible branches in PartialChainSet, and update them by taking actions in GlobalBestActionSet.
But for a PartialChain, which may have 2 actions parallelly to take, how should we do the structure? Do we clone one? whenever we perform the action, we first do the deepcopy.
"""



from copy import deepcopy
from collections import defaultdict
import pandas as pd
import numpy as np
import csv
import json
from pandas import DataFrame
import random
import SeqAttRollOut_RNAfold_only_MFE_4_branch_4_min_Agave_cleared as SeqAttRollOut_RNAfold
import RNA
#import ENTRNA
from util.pseudoknot_free import entrna_main
from util.pseudoknot_free_ori import entrna_main_ori
import itertools





def dp_to_bp(dp):
    a_list = []
    bp_array = np.zeros(len(dp),dtype = int)
    for i in range(len(dp)):
        if dp[i] == "(":
            a_list.append(i)
        if dp[i] == ")":
            bp_array[i] = a_list[-1] + 1
            bp_array[a_list[-1]] = i + 1
            a_list.pop()
    return bp_array


def GenerateLocalPossibleActionSet(PartialChain, RNAOriginalChain):
    """
    For each PartialChain in the PartialChainSet, we generate its feasible actions, and also keep the trace(memo) of the PartialChain itself.
    So we get a LocalPossibleActionSet, which is an attribute of the ParitalChain.
    Also, we append all the actions in LocalPossibleActionSet to GlobalPossibleActionSet.
    """
    LocalPossibleActionSet = SeqAttRollOut_RNAfold.GetPossibleActions(PartialChain, RNAOriginalChain)
    return LocalPossibleActionSet


def UpdateBestWholeChainEverInMemory(BestWholeChainEverInMemory, GlobalPossibleActionSet, OriginalRNAChain, scaler, clf, scaler_ori, clf_ori):
    temp = 0
    append_list = []
    for item in GlobalPossibleActionSet:
        position_of_nucleotide = item.PositionOfNucleotide
        reward = item.rewards
        #print 'rewards', reward, 'vs', BestWholeChainEverInMemory.rewards
        if reward >= BestWholeChainEverInMemory.rewards and item.WholeChainBP != [0 for i in range(0, len(OriginalRNAChain))]:
            #print 'FLAG*****BestWC changed', 'new reward:', reward, 'old reward:', BestWholeChainEverInMemory.rewards
            BestWholeChainEverInMemory.rewards = reward
            BestWholeChainEverInMemory.ChainItself = item.WholeChainBP
            BestWholeChainEverInMemory.Parent = item.PartialChain_parent
            BestWholeChainEverInMemory.ori_WC = item.ori_WC
            temp = temp + 1
        elif item.WholeChainBP == [0 for i in range(0, len(OriginalRNAChain))]:
            pass
        else:
            pass
    if temp == 0:
        print "Sorry, there are no matching partial chain as parent"
        #Here we replace the whole chain and scores of the corresponding next step as fortified's.
        temp_1 =0
        for item in GlobalPossibleActionSet:
            position_of_nucleotide = item.PositionOfNucleotide
            temp_nuc = BestWholeChainEverInMemory.ChainItself[position_of_nucleotide]
            if BestWholeChainEverInMemory.ChainItself[position_of_nucleotide] == 2.5:
                temp_nuc = 2
            elif BestWholeChainEverInMemory.ChainItself[position_of_nucleotide] == 1.5:
                temp_nuc = 1
            #print 'item.BasePair', item.BasePair
            #print 'BestWholeChainEverInMemory.ChainItself[position_of_nucleotide]', BestWholeChainEverInMemory.ChainItself[position_of_nucleotide]
            #print 'temp_nuc', temp_nuc
            #print 'item.WholeChainBP[:(position_of_nucleotide-1)]', item.WholeChainBP[:(position_of_nucleotide-1)]
            #print 'BestWholeChainEverInMemory.ChainItself[:(position_of_nucleotide-1)]', BestWholeChainEverInMemory.ChainItself[:(position_of_nucleotide-1)]
            if item.BasePair == temp_nuc and item.WholeChainBP[:(position_of_nucleotide-1)] == BestWholeChainEverInMemory.ChainItself[:(position_of_nucleotide-1)]:
                temp_1 = temp_1 +1
                new_item = deepcopy(item)
                new_item.WholeChainBP = BestWholeChainEverInMemory.ChainItself
                new_item.ori_WC = BestWholeChainEverInMemory.ChainItself
                #print 'item.Parent', item.PartialChain_parent.ChainItself, 'VVVSSS', 'new_item.PartialChain_parent', new_item.PartialChain_parent.ChainItself
                BestWholeChainEverInMemory.Parent = new_item.PartialChain_parent

                new_rewards_1 = entrna_main(''.join(map(str, OriginalRNAChain)), ''.join(map(str, BP_to_DP(BestWholeChainEverInMemory.ChainItself))), scaler, clf)
                new_rewards_2 = entrna_main_ori(''.join(map(str, OriginalRNAChain)), ''.join(map(str, BP_to_DP(BestWholeChainEverInMemory.ChainItself))), scaler_ori, clf_ori)
                new_rewards_3 = abs(new_rewards_1 - new_rewards_2)

                new_item.rewards = new_rewards_2

                append_list.append(new_item)
        #print 'append_list'
        for itera in append_list:
            #print 'New_action_deepcopy', itera.Nucleotide, itera.BasePair, itera.rewards, itera.PartialChain_parent, itera.WholeChainBP
            GlobalPossibleActionSet.append(itera)
        if temp_1 == 0:
            raise Exception("Sorry, there's no corresponding action to fortified solution.")

    return GlobalPossibleActionSet,BestWholeChainEverInMemory


def RollOut(PartialChain, action, OriginalRNAChain):
    """
    In this function, we input our partial chain, which is the full seq and partial BP chain, and we transfer BP chain into hard constraints.
    """
    PartialChainSim = deepcopy(PartialChain)
    PartialChainSim = SeqAttRollOut_RNAfold.TakeAction(PartialChainSim, action)
    now_position = action.PositionOfNucleotide

    #seq = PartialChainSim.ChainItself
    BPchain = PartialChainSim.BasePairChain
    seq = ''.join(map(str, OriginalRNAChain))

    # Get whole chain with both partial chain constraints and our constraints.
    s = GetWholeStrChain_withConstrainedRNAfold(seq, BPchain)

    #break the str s into pieces, and translate it into BPchain
    WholeChain = PartialChainSim
    new_s = [s[i] for i in range(0, len(s))]
    WholeChain.ori_WC = DP_to_OurBP(new_s)
    new_s = ModifyingWholeChainByRNAfold(new_s, now_position)
    newBP = DP_to_OurBP(new_s)
    WholeChain.BasePairChain = newBP
    return WholeChain


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


def GetWholeStrChain_withConstrainedRNAfold(seq, BPchain):
    """
    In this function, due to iterator problem, we first discuss the length of partial chain.
    And then we get the whole str chain by Constrained RNAfold.
    There are 4 situations we need to consider when adding our constraints to RNAfold in order to disable (...):
    Partial chain ends with ( or (. or (.. or (... ok?
    """
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
    #print 'WC_set', WC_set
    #print 'length of WC_set', len(WC_set)
    removal_list = []
    for item in WC_set:
        #print item[0]
        #print '.'*len(item[0])
        dec = DecideIfAllDot(item[0])
        if dec == True:
            #print 'flag remove'
            #print 'this item' ,item
            removal_list.append(item)
    #print 'removal list:', removal_list
    if removal_list != []:
        for item in removal_list:
            WC_set.remove(item)
    #print 'WC_set updated', WC_set
    #print 'length of updated WC_set', len(WC_set)
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


def DP_to_OurBP(DPChain):
    """
    In this function, we transfer the dot-and-bracket notation into our BasePair notation.
    """
    OurBPNotation = []
    for i in range(len(DPChain)):
        if DPChain[i] == '.':
            OurBPNotation.append(0)
        elif DPChain[i] == '(':
            OurBPNotation.append(1.5)
        elif DPChain[i] == ')':
            OurBPNotation.append(2.5)
    return OurBPNotation


def BP_to_DP(BPChain):
    """
    In this function, we transfer the BasePair notation into dot-and-bracket notation.
    """
    DPNotation = []
    for i in range(len(BPChain)):
        if BPChain[i] == 0:
            DPNotation.append('.')
        elif BPChain[i] == 1.5:
            DPNotation.append('(')
        elif BPChain[i] == 2.5:
            DPNotation.append(')')
        elif BPChain[i] == 1:
            DPNotation.append('x')
        elif BPChain[i] == 2:
            DPNotation.append('y')
    return DPNotation

def UpdateRewardsForAction(PossibleActionSet, OriginalRNAChain, scaler, clf, scaler_ori, clf_ori):
    for action in PossibleActionSet:
        PartialChain = action.PartialChain_parent
        WholeChain = RollOut(PartialChain, action, OriginalRNAChain)
        rewards_foldability_withoutMFE = GetRewards(WholeChain,OriginalRNAChain, scaler, clf)
        rewards_foldability_withMFE = GetRewards_ori(WholeChain,OriginalRNAChain, scaler_ori, clf_ori)
        rewards_diverse = abs(rewards_foldability_withoutMFE - rewards_foldability_withMFE)

        action.rewards = rewards_foldability_withMFE

        action.WholeChainBP = WholeChain.BasePairChain
        action.ori_WC = WholeChain.ori_WC
    return PossibleActionSet


def GenereateTopFourActionSet(GlobalPossibleActionSet):
    """
    We generate the TopActionSet, which is all the top actions with respect to rank in each critics.
    We don't care which action is chosen by which critic.
    """
    temp_GlobalPossibleActionSet = []
    for item in GlobalPossibleActionSet:
        temp_GlobalPossibleActionSet.append(item)
    #print 'temp_GlobalPossibleActionSet', temp_GlobalPossibleActionSet
    TopActionSet_temp = []
    if len(GlobalPossibleActionSet) < 4:
        while temp_GlobalPossibleActionSet != []:
            BestForRewards_top_temp = SelectBestAction(temp_GlobalPossibleActionSet)
            temp_GlobalPossibleActionSet.remove(BestForRewards_top_temp)
            TopActionSet_temp.append(BestForRewards_top_temp)
    else:
        for itera in range(0,4):
            BestForRewards_top_temp = SelectBestAction(temp_GlobalPossibleActionSet)
            temp_GlobalPossibleActionSet.remove(BestForRewards_top_temp)
            TopActionSet_temp.append(BestForRewards_top_temp)
    return TopActionSet_temp


def GetRewards(WholeChain,OriginalRNAChain,scaler, clf):
    """
    In this function we plug in ENTRNA as rewards.
    """
    seq = OriginalRNAChain
    dp = SeqAttRollOut_RNAfold.TransferIntoDP(WholeChain)
    seq_str = ''.join(map(str, seq))
    dp_str = ''.join(map(str, dp))
    #bp = dp_to_bp(dp)
    print('*Get Rewards:*')
    print('seq_str is:',seq_str)
    print('dp_str is:',dp_str)
    foldability = entrna_main(seq_str,dp_str,scaler, clf)
    rewards = foldability
    #print('Now the rewards is',rewards)
    #rewards = random.randrange(0,5)
    return rewards


def GetRewards_ori(WholeChain,OriginalRNAChain,scaler_ori, clf_ori):
    """
    In this function we plug in ENTRNA as rewards.
    """
    seq = OriginalRNAChain
    dp = SeqAttRollOut_RNAfold.TransferIntoDP(WholeChain)
    seq_str = ''.join(map(str, seq))
    dp_str = ''.join(map(str, dp))
    #bp = dp_to_bp(dp)
    print('*Get Rewards Ori:*')
    print('seq_str is:',seq_str)
    print('dp_str is:',dp_str)
    foldability = entrna_main_ori(seq_str,dp_str,scaler_ori, clf_ori)
    rewards = foldability
    #print('Now the rewards is',rewards)
    #rewards = random.randrange(0,5)
    return rewards


def GetRewards1(WholeChain,OriginalRNAChain):
    """
    In this function we plug in ENTRNA as rewards.
    $$$$$ DONT DELETE $$$$$
    """
    seq = OriginalRNAChain
    dp = SeqAttRollOut_RNAfold.TransferIntoDP(WholeChain)
    seq_str = ''.join(map(str, seq))
    dp_str = ''.join(map(str, dp))
    #bp = dp_to_bp(dp)
    print('seq_str is:',seq_str)
    print('dp_str is:',dp_str)
    rewards = entrna_main1(seq_str,dp_str)
    return rewards


def SelectBestAction(PossibleActionSet):
    """
    In this fuction we choose the action with best reward.
    """
    temp = PossibleActionSet[0]
    for item in PossibleActionSet:
        reward = item.rewards
        reward_temp = temp.rewards
        if reward > reward_temp:
            temp = item
    return temp


def WholeChainElementToAction(element):
    if element == 1.5:
        res = 1
    elif element ==2.5:
        res = 2
    else:
        res = element
    return res


def UpdatePartialChainSet(itera, PartialChainSet, RNAOriginalChain, scaler, clf, scaler_ori, clf_ori, BestWholeChainEverInMemory):
    GlobalPossibleActionSet = []
    #print '*initial GPA at each iteration*', GlobalPossibleActionSet
    for pc in PartialChainSet:
        LocalPossibleActionSet = GenerateLocalPossibleActionSet(pc, RNAOriginalChain)
        for lpa in LocalPossibleActionSet:
            #print 'lpa', lpa.Nucleotide, lpa.BasePair, lpa.PositionOfNucleotide
            GlobalPossibleActionSet.append(lpa)
    #print('*GlobalPossibleActionSet:*')
    #for gpa in GlobalPossibleActionSet:
    #    print('Action:', 'Nucleotide:', gpa.Nucleotide, 'BasePair:', gpa.BasePair, 'PositionOfNucleotide:', gpa.PositionOfNucleotide,  'PartialChain_parent_BP_Chain', gpa.PartialChain_parent.BasePairChain)

    UpdateRewardsForAction(GlobalPossibleActionSet, RNAOriginalChain, scaler, clf, scaler_ori, clf_ori)

    #print '*Updated Global Action Set*'
    #for gpa in GlobalPossibleActionSet:
    #    print('Action:', 'Nucleotide:', gpa.Nucleotide, 'PositionOfNucleotide:', gpa.PositionOfNucleotide,'BasePair:', gpa.BasePair ,'AllowableOtherSide:', gpa.AllowableOtherSide, 'Reward', gpa.rewards, 'whole chain', gpa.WholeChainBP, 'ori_WC', gpa.ori_WC)

    #print '*Length of Global Action Set*:', len(GlobalPossibleActionSet)

    GlobalPossibleActionSet,BestWholeChainEverInMemory = UpdateBestWholeChainEverInMemory(BestWholeChainEverInMemory, GlobalPossibleActionSet, RNAOriginalChain, scaler, clf, scaler_ori, clf_ori)


    #print('*GlobalPossibleActionSet_WithFortified:*')
    #for gpa in GlobalPossibleActionSet:
    #    print('Action:', 'Nucleotide:', gpa.Nucleotide, 'PositionOfNucleotide:', gpa.PositionOfNucleotide,'BasePair:', gpa.BasePair ,'AllowableOtherSide:', gpa.AllowableOtherSide, 'Reward', gpa.rewards,  'Parent', gpa.PartialChain_parent.ChainItself, 'WC', gpa.WholeChainBP)

    TopActionSet = GenereateTopFourActionSet(GlobalPossibleActionSet)
    print('*TopActionSet*')
    for ta in TopActionSet:
        print('Action:', 'Nucleotide:', ta.Nucleotide, 'BasePair:', ta.BasePair,  'AllowableOtherSide:', ta.AllowableOtherSide, 'rewards', ta.rewards, 'PartialChain_parent_BP_Chain', ta.PartialChain_parent.BasePairChain)

    PartialChainSet = []
    for action in TopActionSet:
        New_PartialChain= SeqAttRollOut_RNAfold.TakeAction(action.PartialChain_parent, action)
        PartialChainSet.append(New_PartialChain)
    return PartialChainSet, BestWholeChainEverInMemory



def Main_function(OriginalRNAChain,scaler, clf, scaler_ori, clf_ori):
    RNAOriginalChain = OriginalRNAChain
    iteration_length = len(RNAOriginalChain)
    PartialChainSet = []
    BestWholeChainEverInMemory = SeqAttRollOut_RNAfold.BestWholeChainEverInMemory([], 0, 'initial')
    Max_branch_num = 4

    InitialPartialChain = SeqAttRollOut_RNAfold.PartialChain(ChainItself = [], ConnectivityLeft = 1, ConnectivityRight =1, BasePairChain = [], PreservedList = [], MaxAcceptantList = [None for x in range(len(OriginalRNAChain))])
    #print('***InitialPartialChain***',InitialPartialChain.ChainItself, InitialPartialChain.BasePairChain, InitialPartialChain.MaxAcceptantList)
    PartialChainSet.append(InitialPartialChain)

    for itera in range(iteration_length):
        print('**************************Begin Iteration*************************************')
        print('**Iteration**',itera)
        print('**OriginalRNAChain**', OriginalRNAChain)
        PartialChainSet, BestWholeChainEverInMemory = UpdatePartialChainSet(itera, PartialChainSet, RNAOriginalChain, scaler, clf, scaler_ori, clf_ori, BestWholeChainEverInMemory)
        print('***ParitalChainSet:***')
        for item in PartialChainSet:
            print('PartialChain:', 'BP_Chain:', item.BasePairChain, 'DP_Chain', SeqAttRollOut_RNAfold.TransferIntoDP(item))
        print('**************************END OF ITER*****************************************')

    print('**************************Statistics*****************************************')
    print('the length of RNAchain:', len(RNAOriginalChain))
    Our_dp_str_list = []
    print('Foldability:')
    for item in PartialChainSet:
        print('Foldability for this chain:', GetRewards(item, OriginalRNAChain, scaler, clf))
        Our_BP = item.BasePairChain
        Our_BPCon = ''.join(map(str, Our_BP))
        Our_dp = SeqAttRollOut_RNAfold.TransferIntoDP(item)
        Our_dp_str = ''.join(map(str, Our_dp))
        Our_dp_str_list.append(Our_dp_str)
    print('*******************************************************************')
    print('*******************************************************************')

    if len(Our_dp_str_list) == 3:
        Our_dp_str_list.append('NONE')
    elif len(Our_dp_str_list) ==2:
        Our_dp_str_list.append('NONE')
        Our_dp_str_list.append('NONE')
    elif len(Our_dp_str_list) ==1:
        Our_dp_str_list.append('NONE')
        Our_dp_str_list.append('NONE')
        Our_dp_str_list.append('NONE')

    return Our_dp_str_list






