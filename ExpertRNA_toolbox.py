#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is the toolbox to save functions to call.
"""

from copy import deepcopy
from collections import defaultdict
import pandas as pd
import numpy as np
import csv
import json
from pandas import DataFrame
import random
from re import finditer

#-------------------------Classes definition-------------------------
class Action():
    def __init__(self, Nucleotide, BasePair, PositionOfNucleotide, PartialChain_parent):
        """
        For each action, we have four elements of it:
            (1) The nucleotide itself
            (2) Base pair:
                0: Unpair.
                1: Open a base pair.
                2: Close a base pair.
            (3) PositionOfNucleotide defines at which position we pin the nucleotide on the given original chain.
	        (4) Rewards for each critics: the reward is in the form of {'ExpertName': [reward, folder]}, and each expert we pick the best performing WC from each folder
	        (5) PartialChain_parent: The parent branch for this action.
	        (6) ori_WC: wc after rollout, in the form of {'FolderName': WC}

        """
        self.Nucleotide = Nucleotide  #change to LabelOfNucleotide
        self.BasePair = BasePair  #change to 
        self.PositionOfNucleotide = PositionOfNucleotide
        self.rewards = {}
        self.PartialChain_parent = PartialChain_parent
        self.ori_WC = {}



class Solution():  #change it to solution
    def __init__(self, ChainItself, BasePairChain):
        """
        This 'Solution' used to be called PartialChain.
        For each action, we have four elements of it:
            (1) The partial chain itself
            (2) Store the BP information
            (4) ori_WC: it is the fortified wholechain saved
                it is in the form of a marix, row is by the num of folders and col is by the experts
        """
        self.ChainItself = ChainItself
        self.BasePairChain = BasePairChain
        #self.MaxAcceptantList = MaxAcceptantList
        self.ori_WC = 'NNN'


#---------------------------------------------------------------------------------------------------------------------

# returns a vector containing each nucleotide's pair
def parse_dot_bracket(input):
    output = np.full(len(input), -1)
    # I'm not sure this is the most efficent way to do this, but I'm lazy.
    more = True
    while more:
        more = False

        # finds matched parenthesis
        for x in finditer(r"\([^()]*\)", input):
            more = True
            output[x.start()] = x.end() - 1
            output[x.end() - 1] = x.start()

            input = input[0:x.start()] + "." + input[x.start() + 1:x.end() - 1] + "." + input[x.end():]

    return output

def TurnWholeChainToPCPairingInfo(WholeChain, PartialChain):
    #In this function, we grab the pairing info of partial chain from the whole chain.
    #The WC can be WC saved in PC from last step

    print('TurnWholeChainToPCPairingInfo:')
    print('WholeChain:', WholeChain)
    print('PartialChain:', PartialChain)

    AllOpenPosInPC = []
    for i in range(0, len(PartialChain)):
        if PartialChain[i] == 1 or PartialChain[i] == 1.5:
            AllOpenPosInPC.append(i)
    print('AllOpenPosInPC:', AllOpenPosInPC)

    # Turn [1,1.5,0,0,0,0,2.5,...] form into [[1,6]] form
    pair_table = parse_dot_bracket(WholeChain).tolist()
    pair_table_points = []
    for i in range(0, len(pair_table)):
        if pair_table[i] != -1 and i < pair_table[i]:
            pair_table_points.append([i, pair_table[i]])
    print('pair_table_points:', pair_table_points)

    # grab pairing info only for opens in pc
    # the form  is  [[1,6],[2,5],...] like so
    PCPairingInfo = []
    for p in AllOpenPosInPC:
        for item in pair_table_points:
            if item[0] == p:
                PCPairingInfo.append(item)

    return PCPairingInfo

def dp_to_bp(dp):
    a_list = []
    bp_array = np.zeros(len(dp), dtype=int)
    for i in range(len(dp)):
        if dp[i] == "(":
            a_list.append(i)
        if dp[i] == ")":
            bp_array[i] = a_list[-1] + 1
            bp_array[a_list[-1]] = i + 1
            a_list.pop()
    return bp_array

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



def DecideIfAllDot(s):
    temp = 0
    for i in range(0, len(s)):
        if s[i] != '.':
            temp = temp + 1
        else:
            pass
    if temp == 0:
        decision = True
    else:
        decision = False
    return decision



#---------------------------------------------------------------------------------------------------------------------

def UpdateBPChain(BPChain, ActionChosen):
    """
    This function updates the BP chain state after we attach a new one, it transfers uncomplete to newly complete: 1 to 1.5, 2 to 2.5
    """
    print('--*--UpdateBPChain--*--')
    # print('test_ori_BP',BPChain)
    # NewBPChain = deepcopy(BPChain)
    NewBPChain = BPChain
    if ActionChosen.BasePair == 0:
        pass
    elif ActionChosen.BasePair == 1:
        pass
    elif ActionChosen.BasePair == 2:
        # print('FLAG ActionUpdate as 2')
        # print('Length of BPChain:',len(NewBPChain))
        # print('testNewBPChainOriginal:',NewBPChain)
        PositionOfAttachedNucleotide = ActionChosen.PositionOfNucleotide
        PositionOfLastOpenInPartialChain = len(NewBPChain) - 1 - NewBPChain[::-1].index(1)
        print('PositionOfAttachedNucleotide', PositionOfAttachedNucleotide)
        print('PositionOfLastOpenInPartialChain', PositionOfLastOpenInPartialChain)
        print('length of NewBPChain', len(NewBPChain))
        NewBPChain[PositionOfAttachedNucleotide] = 2.5
        NewBPChain[PositionOfLastOpenInPartialChain] = 1.5
    # print('testNewBPChainUpdated:',NewBPChain)
    return NewBPChain


def TakeAction(PartialChain, ActionChosen):
    """
    Here we take the action chosen and update the partial chain.
    $$$We also need to update BPchain.
    """
    print('--*--TakeAction--*--')
    NewPartialChain = deepcopy(PartialChain)
    # NewPartialChain = PartialChain
    NewPartialChain.ChainItself.append(ActionChosen.Nucleotide)
    NewPartialChain.BasePairChain.append(ActionChosen.BasePair)
    # print NewPartialChain.ChainItself
    # print NewPartialChain.BasePairChain
    NewPartialChain.BasePairChain = UpdateBPChain(NewPartialChain.BasePairChain, ActionChosen)
    NewPartialChain.ori_WC = ActionChosen.rewards
    # print('NewPartialChain_Nucleotide', NewPartialChain.ChainItself, 'NewPartialChain_BasePairChain',NewPartialChain.BasePairChain)
    return NewPartialChain

# ---------------------------------------------------------------------------------------------------------------------

def GetPossibleActions(Solution, RNAOriginalChain, min_dbp):
    """
    This function is used to generate feasible action for the solution/PartialChain now.
    """
    PossibleActionSet = []
    piece = Solution.ChainItself
    BPpiece = Solution.BasePairChain
    #MaxAcceptantList = Solution.MaxAcceptantList
    p = len(piece) #position of this nucleotide
    #this_nucleotide = RNAOriginalChain[p]
    #min_dbp = min_dbp   #min distance btw base pair open and close

    # Attach possible actions
    NucleotideChoiceList = []
    #AttachToDirectionChoiceList = []
    #AllowableOtherSideChoiceList = []
    PositionOfNucleotideChoiceList = []

    NucleotideChoiceList.append(RNAOriginalChain[p])
    #AttachToDirectionChoiceList.append('left')
    #AllowableOtherSideChoiceList.append(1)
    PositionOfNucleotideChoiceList.append(p)

    # About base pair info
    """
    In the base pair chain, we define 
        uncomplete open as 1, complete open as 1.5,
        uncomplete close as 2, complete close as 2.5.
            
    """
    CountOfOpenInPartialChainNotComplete = BPpiece.count(1)
    CountOfCloseInPartialChainNotComplete = BPpiece.count(2)
    if CountOfCloseInPartialChainNotComplete != 0:
        raise Exception("Sorry, CountOfCloseInPartialChainNotComplete != 0")
    # print('CountOfOpenInPartialChainNotComplete',CountOfOpenInPartialChainNotComplete,'CountOfCloseInPartialChainNotComplete',CountOfCloseInPartialChainNotComplete)

    if CountOfOpenInPartialChainNotComplete != 0:
        # print('Situation1')
        PositionOfFirstOpenInPartialChain = int(BPpiece.index(1))
        #PositionOfFirstCloseInPartialChain = int(BPpiece.index(2))
        PositionOfLastOpenInPartialChain = len(BPpiece) - BPpiece[::-1].index(1) - 1
        #PositionOfLastCloseInPartialChain = len(BPpiece) - BPpiece[::-1].index(2) - 1
    elif CountOfOpenInPartialChainNotComplete == 0:
        # print('Situation2')
        PositionOfFirstOpenInPartialChain = 'non'
        #PositionOfFirstCloseInPartialChain = int(BPpiece.index(2))
        PositionOfLastOpenInPartialChain = 'non'
        #PositionOfLastCloseInPartialChain = int(len(BPpiece) - BPpiece[::-1].index(2) - 1)
    else:
        pass

    #begin
    BasePairChoiceList = []
    #Can we attach unpair?
    if PositionOfLastOpenInPartialChain != 'non':
        # and UpperboundForInterval-PositionOfNucleotideChoiceList[i] >= CountOfOpenInPartialChainNotCompleteWithinInterval:#
        if CheckIfThisOneIsUnpairWillTheRemaingingFeasible(RNAOriginalChain, PositionOfLastOpenInPartialChain, p, BPpiece, min_dbp) == True:
            BasePairChoiceList.append(0)
    elif PositionOfLastOpenInPartialChain == 'non':
        BasePairChoiceList.append(0)
    #Can we attach close?
    if PositionOfLastOpenInPartialChain == 'non':
        pass
    else:
        if CheckIfNucleotidesCanBePaired(PositionOfLastOpenInPartialChain, p, RNAOriginalChain) ==1:
            BasePairChoiceList.append(2)
    #Can we attach open?
    if p + min_dbp + 1 <= len(RNAOriginalChain):
        print('flag1')
        if CheckIfThisOneIOpenWillTheRemaingingFeasible(RNAOriginalChain, PositionOfLastOpenInPartialChain, p, BPpiece, min_dbp) == True:
            print('flag2')
            BasePairChoiceList.append(1)
    
    for j in range(len(BasePairChoiceList)):
        PossibleActionSet.append(Action(Nucleotide=RNAOriginalChain[p], BasePair=BasePairChoiceList[j],PositionOfNucleotide=p, PartialChain_parent = Solution))

    return PossibleActionSet



def CheckIfNucleotidesCanBePaired(PositionOfOpenIndex, PositionOfCloseIndex, NucleotideChain):
    #print('open:', PositionOfOpenIndex, 'close:', PositionOfCloseIndex)
    PositionOfOpen = NucleotideChain[PositionOfOpenIndex]
    #print('PositionOfOpen is:',PositionOfOpen)
    PositionOfClose = NucleotideChain[PositionOfCloseIndex]
    #print('PositionOfClose is:',PositionOfClose)
    if PositionOfOpen == 'A' and PositionOfClose == 'U':
        check = 1
    elif PositionOfOpen == 'G':
        if PositionOfClose == 'C' or PositionOfClose == 'U':
            check = 1
        else:
            check = 0
    elif PositionOfOpen == 'U':
        if PositionOfClose == 'A' or PositionOfClose == 'G':
            check = 1
        else:
            check = 0
    elif PositionOfOpen == 'C' and PositionOfClose == 'G':
        check = 1
    else:
        check = 0
    return check


def CheckIfThisOneIOpenWillTheRemaingingFeasible(RNAOriginalChain, PositionOfLastOpenInPartialChain, PositionOfThisNucleotide, BPpiece, min_dbp):
    """
    ***NEW***
    In here, we not only consider about how many  potential acceptant to come,but also considering all the open which may consume acceptants.
    Same thought, if we added this '(', can we find a feasible solution for the remaining?
    """
    RNAOriginalChain_copy = deepcopy(RNAOriginalChain)
    #RNAOriginalChain_copy = RNAOriginalChain_copy.split(",")
    BPpiece_copy = deepcopy(BPpiece)
    # print('*The BPpiece Check:', BPpiece_copy)
    BPpiece_copy.append(1)
    # print('*The BPpiece Check Updated:', BPpiece_copy)
    #print('RNAOriginalChain_copy:', RNAOriginalChain_copy)
    #print('BPpiece_copy:', BPpiece_copy)
    checkDic = {}
    for i in range(len(BPpiece_copy) - 1, -1, -1):
        if BPpiece_copy[i] == 1:
            checkDic[str(i)] = 'not yet paired'
            NucleotideType = RNAOriginalChain_copy[i]
            if NucleotideType == 'A':
                AcceptantType = ['U']
            elif NucleotideType == 'G':
                AcceptantType = ['C', 'U']
            elif NucleotideType == 'U':
                AcceptantType = ['A', 'G']
            elif NucleotideType == 'C':
                AcceptantType = ['G']
            if 'Taken' in RNAOriginalChain_copy:
                LastTakenPosition = len(RNAOriginalChain_copy) - RNAOriginalChain_copy[::-1].index('Taken') - 1
            else:
                LastTakenPosition = 0
            # print('LastTakenPosition',LastTakenPosition)
            for j in range(max(PositionOfThisNucleotide + min_dbp + 1, LastTakenPosition + 1), len(RNAOriginalChain_copy)):
                for Acceptant in AcceptantType:
                    if RNAOriginalChain[j] == Acceptant:
                        checkDic[str(i)] = 'paired'
                        RNAOriginalChain_copy[i] = 'TakenOpen'
                        RNAOriginalChain_copy[j] = 'Taken'
                        break
                else:
                    continue
                break
        else:
            pass
            # print('*RNAchain now_temp:',RNAOriginalChain_copy)
    if checkDic == {}:
        check = True
    else:
        print('checkDic:', checkDic)
        check = all(value == 'paired' for value in checkDic.values())
    # print('*checkDic:',checkDic)
    # print('*RNAchain now:',RNAOriginalChain_copy)
    # print('*check is:',check)
    return check


def CheckForInterval(PositionOfThisNucleotide, RNAOriginalChain, MaxAcceptantList, PositionOfLastOpenInPartialChain):
    """
    In this function, we look for upperbound for the smallest interval surrounding our position now.
    The lower bound is just PositionOfLastOpenInPartialChain.
    NOT IN USE.
    """
    if PositionOfLastOpenInPartialChain == 'non':
        Upperbound = len(RNAOriginalChain)
    elif MaxAcceptantList[PositionOfLastOpenInPartialChain] == None:
        Upperbound = len(RNAOriginalChain)
    else:
        Upperbound = MaxAcceptantList[PositionOfLastOpenInPartialChain]
    return Upperbound


def CheckIfThisOneIsUnpairWillTheRemaingingFeasible(RNAOriginalChain, PositionOfLastOpenInPartialChain, PositionOfThisNucleotide, BPpiece, min_dbp):
    """
    ***NEW***
    In this function we are trying to find one feasible solution, under the assumption that if this position is unpair.
    """
    RNAOriginalChain_copy = deepcopy(RNAOriginalChain)
    #RNAOriginalChain_copy = RNAOriginalChain_copy.split(",")
    BPpiece_copy = deepcopy(BPpiece)
    RNAOriginalChain_copy[PositionOfThisNucleotide] = 'Unpaired'
    # print('The BPpiece Check:', BPpiece_copy)
    checkDic = {}
    LastTakenPosition = 0
    for i in range(len(BPpiece_copy) - 1, -1, -1):
        if BPpiece_copy[i] == 1:
            checkDic[str(i)] = 'not yet paired'
            NucleotideType = RNAOriginalChain_copy[i]
            if NucleotideType == 'A':
                AcceptantType = ['U']
            elif NucleotideType == 'G':
                AcceptantType = ['C', 'U']
            elif NucleotideType == 'U':
                AcceptantType = ['A', 'G']
            elif NucleotideType == 'C':
                AcceptantType = ['G']
            if 'Taken' in RNAOriginalChain_copy:
                LastTakenPosition = len(RNAOriginalChain_copy) - RNAOriginalChain_copy[::-1].index('Taken') - 1
            else:
                LastTakenPosition = 0
            # print('LastTakenPosition',LastTakenPosition)
            for j in range(max(PositionOfLastOpenInPartialChain + 1 + min_dbp, LastTakenPosition + 1, PositionOfThisNucleotide + 1),len(RNAOriginalChain_copy)):
                for Acceptant in AcceptantType:
                    if RNAOriginalChain[j] == Acceptant:
                        checkDic[str(i)] = 'paired'
                        RNAOriginalChain_copy[i] = 'TakenOpen'
                        RNAOriginalChain_copy[j] = 'Taken'
                        break
                else:
                    continue
                break
            # print('RNAchain now',RNAOriginalChain_copy)
        else:
            pass

    if checkDic == {}:
        check = True
    else:
        check = all(value == 'paired' for value in checkDic.values())
    # print('checkDic:',checkDic)
    # print('RNAchain now:',RNAOriginalChain_copy)
    # print('check is:',check)
    return check

