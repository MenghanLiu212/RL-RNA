"""
Implementation of Alg#2 of RL-RNA
We are trying to attach possible actions to the set possible action set.

$$$$$(Oct.6 15:13)
Add dot and bracket notation to the basepair info.


$$$$$(Oct.8 11:00)
Add constriants for base pair choices.

$$$$$$(Oct.20)
Change the nucleotide attaching order as 5' to 3', now we only consider the BP information.
"""

from copy import deepcopy
from collections import defaultdict
import pandas as pd
import numpy
import csv
import json
from pandas import DataFrame
import random
import Main_SeqAttRollOut_RNAfold_only_MFENFE_2X2_branches_4_min_Agave_5s as Our_alg_Main





class Action():
    def __init__(self, Nucleotide, AttachToDirection, AllowableOtherSide, BasePair, PositionOfNucleotide, rewards_1, rewards_2, WholeChainBP, PartialChain_parent):
        """
        For each action, we have four elements of it:
            (1) The nucleotide itself
            (2) Attach direction to the existing partial chain. It's left or right.
            (3) If the other direction is allowable, if not, it would be the 5' end or 3' end. o as no, 1 as yes.
            (4) Base pair:
                0: No base pair.
                1: Open a base pair.
                2: Close a base pair.
            (5) PositionOfNucleotide defines at which position we pin the nucleotide on the given original chain.
	    (6) Rewards for each critics
	    (7) WholeChainBP: whole chain after rollout
	    (8) PartialChain_parent: The parent branch for this action.
        """
        self.Nucleotide = Nucleotide
        self.AttachToDirection = AttachToDirection
        self.AllowableOtherSide = AllowableOtherSide
        self.BasePair = BasePair
        self.PositionOfNucleotide = PositionOfNucleotide
        self.rewards_1 = rewards_1
        self.rewards_2 = rewards_2
        self.WholeChainBP = WholeChainBP
        self.PartialChain_parent = PartialChain_parent
        self.ori_WC = 'Not yet decided'

class PartialChain():
    def __init__(self, ChainItself, ConnectivityLeft, ConnectivityRight, BasePairChain, PreservedList, MaxAcceptantList):
        """
        For each action, we have four elements of it:
            (1) The partial chain itself
            (2) Whether the left side of the partial chain is allowed to connect other nicleotide.
            (3) Whether the right side of the partial chain is allowed to connect other nicleotide.
            (4) Store the BP information
            (5) PreservedList: This is for action preservation, we store each element as a length4 list as [cause,BP[cause],result,BP[result]]
        """
        self.ChainItself = ChainItself
        self.ConnectivityLeft = ConnectivityLeft
        self.ConnectivityRight = ConnectivityRight
        self.BasePairChain = BasePairChain
        self.PreservedList = PreservedList
        self.MaxAcceptantList = MaxAcceptantList
        self.ori_WC = 'Not yet decided'

class BestWholeChainEverInMemory():
    def __init__(self, ChainItself, metric, rewards, Parent):
        """
        This is storing best WC ever in memory for each Metric(expert).
        """
        self.ChainItself = ChainItself
        self.metric = metric
        self.rewards = rewards
        self.Parent = Parent
        self.ori_WC = 'Not yet decided'
#---------------------------------------------------------------------------------------------------------------------

def GetPossibleActions(PartialChain, RNAOriginalChain):
    PossibleActionSet = []
    Positions = []
    ori_chain = RNAOriginalChain
    piece = PartialChain.ChainItself
    BPpiece = PartialChain.BasePairChain
    PreservedList = PartialChain.PreservedList
    MaxAcceptantList =PartialChain.MaxAcceptantList
    p = len(piece)
    
    #Attach possible actions
    NucleotideChoiceList = []
    AttachToDirectionChoiceList = []
    AllowableOtherSideChoiceList = []
    PositionOfNucleotideChoiceList = []
        
    NucleotideChoiceList.append(RNAOriginalChain[p])
    AttachToDirectionChoiceList.append('left')
    AllowableOtherSideChoiceList.append(1)
    PositionOfNucleotideChoiceList.append(p)
               
    #About base pair info
    """
    In the base pair chain, we define 
        uncomplete open as 1, complete open as 1.5,
        uncomplete close as 2, complete close as 2.5.
            
    """   
    CountOfOpenInPartialChainNotComplete = BPpiece.count(1)
    CountOfCloseInPartialChainNotComplete = BPpiece.count(2)
    #print('CountOfOpenInPartialChainNotComplete',CountOfOpenInPartialChainNotComplete,'CountOfCloseInPartialChainNotComplete',CountOfCloseInPartialChainNotComplete)
    LeftPositionOfPartialChain = 0
    RightPositionOfPartialChain = p
    
    if CountOfOpenInPartialChainNotComplete!=0 and CountOfCloseInPartialChainNotComplete!= 0:
        #print('Situation1')
        PositionOfFirstOpenInPartialChain = int(BPpiece.index(1))
        PositionOfFirstCloseInPartialChain = int(BPpiece.index(2))
        PositionOfLastOpenInPartialChain = len(BPpiece) - BPpiece[::-1].index(1) -1
        PositionOfLastCloseInPartialChain = len(BPpiece) - BPpiece[::-1].index(2) -1
        #NumberOfRemainingPositionsRightToPartialChain = len(BPpiece)-RightPositionOfPartialChain
        #NumberOfRemainingPositionsLeftToPartialChain = LeftPositionOfPartialChain-1
    elif CountOfOpenInPartialChainNotComplete==0 and CountOfCloseInPartialChainNotComplete!= 0:
        #print('Situation2')
        PositionOfFirstOpenInPartialChain = 'non'
        PositionOfFirstCloseInPartialChain = int(BPpiece.index(2))
        PositionOfLastOpenInPartialChain = 'non'
        PositionOfLastCloseInPartialChain = int(len(BPpiece) - BPpiece[::-1].index(2) -1)
    elif CountOfOpenInPartialChainNotComplete!=0 and CountOfCloseInPartialChainNotComplete== 0:
        #print('Situation3')
        PositionOfFirstOpenInPartialChain = BPpiece.index(1)
        PositionOfFirstCloseInPartialChain ='non'
        PositionOfLastOpenInPartialChain = int(len(BPpiece) - BPpiece[::-1].index(1) -1)
        PositionOfLastCloseInPartialChain = 'non'
    else:
        #print('Situation4')
        PositionOfFirstOpenInPartialChain = 'non'#
        PositionOfFirstCloseInPartialChain = 'non'
        PositionOfLastOpenInPartialChain = 'non'
        PositionOfLastCloseInPartialChain = 'non'#
       #####print('PositionOfFirstOpenInPartialChain',PositionOfFirstOpenInPartialChain,'PositionOfFirstCloseInPartialChain',PositionOfFirstCloseInPartialChain,'PositionOfLastOpenInPartialChain',PositionOfLastOpenInPartialChain,'PositionOfLastCloseInPartialChain',PositionOfLastCloseInPartialChain)
    if 1==1: #just for identing
        i =len(NucleotideChoiceList)-1
        print 'i', i
        #####print('PositionOfNucleotideChoiceList[i]',PositionOfNucleotideChoiceList[i])
        BasePairChoiceList = []
        if PositionOfNucleotideChoiceList[i] >= RightPositionOfPartialChain: #attach to the right part of partial chain
            ###Interval
            UpperboundForInterval = CheckForInterval(PositionOfNucleotideChoiceList[i], RNAOriginalChain, MaxAcceptantList, PositionOfLastOpenInPartialChain)
            if PositionOfLastOpenInPartialChain != 'non':
                #pieceWithinInterval = piece[PositionOfLastOpenInPartialChain:UpperboundForInterval]
                BPpieceWithinInterval = BPpiece[PositionOfLastOpenInPartialChain:UpperboundForInterval]
            else:
                #pieceWithinInterval = piece[0:UpperboundForInterval]
                BPpieceWithinInterval = BPpiece[0:UpperboundForInterval]
            CountOfOpenInPartialChainNotCompleteWithinInterval = BPpieceWithinInterval.count(1)
            #print('UpperboundForInterval',UpperboundForInterval)
            ###Not pairing choice
            if PositionOfLastOpenInPartialChain!='non': 
            #and UpperboundForInterval-PositionOfNucleotideChoiceList[i] >= CountOfOpenInPartialChainNotCompleteWithinInterval:#
		if CheckIfThisOneIsUnpairWillTheRemaingingFeasible(PartialChain, ori_chain, PositionOfLastOpenInPartialChain, PositionOfNucleotideChoiceList[i], BPpiece) == True:
                    BasePairChoiceList.append(0)
            elif PositionOfLastOpenInPartialChain=='non':
                BasePairChoiceList.append(0)
            ###Condtions
            if PositionOfLastOpenInPartialChain == PositionOfLastCloseInPartialChain:
                ConditionAllNon = 1
            else:
                ConditionAllNon = 0
            if ConditionAllNon == 0:
                if PositionOfLastOpenInPartialChain!='non' and PositionOfLastCloseInPartialChain!='non':  
                    if PositionOfLastOpenInPartialChain > PositionOfLastCloseInPartialChain:
                        Condition1 = 1
                    else:
                        Condition1=0
                else:
                    Condition1=0
                if PositionOfLastOpenInPartialChain!='non' and PositionOfLastCloseInPartialChain=='non': 
                    Condition1_1 = 1
                else:
                    Condition1_1 = 0
                if PositionOfLastOpenInPartialChain!='non' and PositionOfLastCloseInPartialChain!='non': 
                    if PositionOfLastOpenInPartialChain <= PositionOfLastCloseInPartialChain:
                        Condition2 = 1
                    else:
                        Condition2=0
                else:
                    Condition2=0
                if PositionOfLastOpenInPartialChain=='non' and PositionOfLastCloseInPartialChain!='non': 
                    Condition2_1 = 1
                else:
                    Condition2_1 = 0
            ###Start for deciding
            if ConditionAllNon == 0:
                if (Condition1) or (Condition1_1):
                #print('flag_cond1*****')
                    if PositionOfNucleotideChoiceList[i]-PositionOfLastOpenInPartialChain > 4:
                        #print('PositionOfNucleotideChoiceList[i]',PositionOfNucleotideChoiceList[i])
                        if CheckIfNucleotidesCanBePaired(PositionOfLastOpenInPartialChain, PositionOfNucleotideChoiceList[i],RNAOriginalChain) ==1:
                            BasePairChoiceList.append(2)
                    if len(ori_chain)-PositionOfNucleotideChoiceList[i]-1 > 4+CountOfOpenInPartialChainNotComplete:
                        #print('flag*****')
                        #tempList = CalculateNumOfPotentialAcceptant(PositionOfNucleotideChoiceList[i], piece, BPpiece, RNAOriginalChain, MaxAcceptantList, UpperboundForInterval)
                        #MaxAcceptant = tempList[1]
                        #PartialChain.MaxAcceptantList[PositionOfNucleotideChoiceList[i]] = MaxAcceptant
                        #if tempList[0] > 0:
                        #    BasePairChoiceList.append(1)
                        #use another method for deciding open
                        if CheckIfThisOneIOpenWillTheRemaingingFeasible(PartialChain,ori_chain,PositionOfLastOpenInPartialChain,PositionOfNucleotideChoiceList[i], BPpiece) == True:
                            BasePairChoiceList.append(1)
                elif (Condition2) or (Condition2_1):
                    if PositionOfFirstCloseInPartialChain > 4+CountOfCloseInPartialChainNotComplete:
                        if CountOfOpenInPartialChainNotComplete > 0:
                            BasePairChoiceList.append(2)
                    if len(ori_chain)-PositionOfNucleotideChoiceList[i] > 4+1:
                        #if RNAOriginalChain[PositionOfNucleotideChoiceList[i]] != 'C':
                        tempList = CalculateNumOfPotentialAcceptant(PositionOfNucleotideChoiceList[i], piece, BPpiece, RNAOriginalChain, MaxAcceptantList, UpperboundForInterval)
                        MaxAcceptant = tempList[1]
                        PartialChain.MaxAcceptantList[PositionOfNucleotideChoiceList[i]] = MaxAcceptant
                        if tempList[0] > 0:
                            BasePairChoiceList.append(1)
            elif ConditionAllNon == 1:
                #####print('CountOfOpenInPartialChainNotComplete:',CountOfOpenInPartialChainNotComplete)
                #####print('CountOfCloseInPartialChainNotComplete',CountOfCloseInPartialChainNotComplete)
                if PositionOfNucleotideChoiceList[i] > 1+4:
                    if CountOfOpenInPartialChainNotComplete > 0:
                            BasePairChoiceList.append(2)
                if len(ori_chain)-PositionOfNucleotideChoiceList[i] > 1+4:
                    #if RNAOriginalChain[PositionOfNucleotideChoiceList[i]] != 'C':
		    #print('flagAllNon*****')
                    tempList = CalculateNumOfPotentialAcceptant(PositionOfNucleotideChoiceList[i], piece, BPpiece, RNAOriginalChain, MaxAcceptantList, UpperboundForInterval)
                    MaxAcceptant = tempList[1]
                    PartialChain.MaxAcceptantList[PositionOfNucleotideChoiceList[i]] = MaxAcceptant
                    #print('tempList[0]',tempList[0])
                    if tempList[0] > 0:
                        BasePairChoiceList.append(1)       
    #####print('The length of BasePairChoiceList', len(BasePairChoiceList))
    #Possible_temp = []             
    for j in range(len(BasePairChoiceList)):
        PossibleActionSet.append(Action(Nucleotide=NucleotideChoiceList[i], AttachToDirection=AttachToDirectionChoiceList[i], AllowableOtherSide = AllowableOtherSideChoiceList[i], BasePair = BasePairChoiceList[j], PositionOfNucleotide=PositionOfNucleotideChoiceList[i], rewards_1 = 'not deicided yet', rewards_2 = 'not deicided yet', WholeChainBP = 'WC not deicided yet', PartialChain_parent = PartialChain))
        """
	Possible_temp.append(BasePairChoiceList[j])
    if BestWholeChainEverInMemory != []:
	temp = BestWholeChainEverInMemory[PositionOfNucleotideChoiceList[i]]
	print('temp now?',temp)
	if temp == 1.5 or temp==1:
	    temp = 1
	elif temp == 2.5 or temp ==2 :
	    temp = 2
	elif temp == 0:
	    temp = 0
	#decide if the best is feasible
	if temp in Possible_temp:
	    PossibleActionSet.append(Action(Nucleotide=NucleotideChoiceList[i], AttachToDirection=AttachToDirectionChoiceList[i], AllowableOtherSide = 'xxxIMTHEBESTxxx', BasePair = temp, PositionOfNucleotide=PositionOfNucleotideChoiceList[i], rewards = BestWholeChainEverInMemoryScore, WholeChainBP = BestWholeChainEverInMemory))
        """
    return PossibleActionSet


def CheckIfNucleotidesCanBePaired(PositionOfOpenIndex, PositionOfCloseIndex, NucleotideChain):
    PositionOfOpen=NucleotideChain[PositionOfOpenIndex]
    #print('PositionOfOpen is:',PositionOfOpen)
    PositionOfClose=NucleotideChain[PositionOfCloseIndex]
    #print('PositionOfClose is:',PositionOfClose)
    if PositionOfOpen == 'A' and PositionOfClose == 'U':
        check=1
    elif PositionOfOpen == 'G':
        if PositionOfClose == 'C' or PositionOfClose == 'U':
            check=1
        else:
            check=0
    elif PositionOfOpen == 'U':
        if PositionOfClose == 'A' or PositionOfClose == 'G':
            check=1
        else:
            check=0
    elif PositionOfOpen == 'C' and PositionOfClose == 'G':
        check=1
    else:
        check=0
    return check 


def CheckIfThisOneIOpenWillTheRemaingingFeasible(PartialChain, RNAOriginalChain, PositionOfLastOpenInPartialChain, PositionOfThisNucleotide, BPpiece):
    """
    ***NEW***
    In here, we not only consider about how many  potential acceptant to come,but also considering all the open which may consume acceptants.
    Same thought, if we added this '(', can we find a feasible solution for the remaining?
    """
    RNAOriginalChain_copy = deepcopy(RNAOriginalChain)
    BPpiece_copy = deepcopy(BPpiece)
    #print('*The BPpiece Check:', BPpiece_copy)
    BPpiece_copy.append(1)
    #print('*The BPpiece Check Updated:', BPpiece_copy)
    checkDic = {}
    for i in range(len(BPpiece_copy)-1,-1,-1):
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
            if 'Taken'in RNAOriginalChain_copy:
                LastTakenPosition = len(RNAOriginalChain_copy) - RNAOriginalChain_copy[::-1].index('Taken') -1    
            else:
                LastTakenPosition = 0
            #print('LastTakenPosition',LastTakenPosition)
            for j in range(max(PositionOfThisNucleotide+4+1,LastTakenPosition+1), len(RNAOriginalChain_copy)):
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
            #print('*RNAchain now_temp:',RNAOriginalChain_copy)
    if checkDic == {}:
        check = True
    else:
        check = all(value == 'paired' for value in checkDic.values())
    #print('*checkDic:',checkDic)
    #print('*RNAchain now:',RNAOriginalChain_copy)
    #print('*check is:',check)
    return check


def CalculateNumOfPotentialAcceptant(PositionOfThisNucleotide, piece, BPpiece, RNAOriginalChain, MaxAcceptantList, UpperboundForInterval):
    """
    NOTICE: Not using now.
    """
    NucleotideType = RNAOriginalChain[PositionOfThisNucleotide]
    if NucleotideType == 'A':
        AcceptantType = ['U']
    elif NucleotideType == 'G':
        AcceptantType = ['C', 'U']
    elif NucleotideType == 'U':
        AcceptantType = ['A', 'G']
    elif NucleotideType == 'C':
        AcceptantType = ['G']
    temp1 = 0
    for i in range(len(piece)):
        if piece[i] == NucleotideType and BPpiece[i]== 1 :
            temp1 = temp1 + 1
    NumOfOpenAsSameNucleotideType = temp1
    temp2 = 0
    tempMaxList = []
    for item in AcceptantType:
        for i in range(5+PositionOfThisNucleotide, UpperboundForInterval):#
            if RNAOriginalChain[i] == item:
                temp2= temp2 + 1
                tempMaxList.append(i)
    NumOfCloseAsCorrespondingNucleotideType = temp2
    if tempMaxList != []:
        MaxAcceptant = max(tempMaxList)
    else:
        MaxAcceptant = None
    #NumOfPotentialAcceptant = NumOfCloseAsCorrespondingNucleotideType - NumOfOpenAsSameNucleotideType
    if PositionOfThisNucleotide > 0 and MaxAcceptantList[PositionOfThisNucleotide-1]!=None:
        NumOfPotentialAcceptant = MaxAcceptantList[PositionOfThisNucleotide-1] - NumOfOpenAsSameNucleotideType
    elif PositionOfThisNucleotide == 0:
        NumOfPotentialAcceptant = NumOfCloseAsCorrespondingNucleotideType
    elif MaxAcceptantList[PositionOfThisNucleotide-1]==None:
        NumOfPotentialAcceptant = NumOfCloseAsCorrespondingNucleotideType
    #print('NumOfPotentialAcceptant',NumOfPotentialAcceptant)
    return [NumOfPotentialAcceptant,MaxAcceptant]


def CheckForInterval(PositionOfThisNucleotide, RNAOriginalChain, MaxAcceptantList, PositionOfLastOpenInPartialChain):
    """
    In this function, we look for upperbound for the smallest interval surrounding our position now.
    The lower bound is just PositionOfLastOpenInPartialChain.
    """
    if PositionOfLastOpenInPartialChain == 'non':
        Upperbound = len(RNAOriginalChain)
    elif MaxAcceptantList[PositionOfLastOpenInPartialChain] == None:
        Upperbound = len(RNAOriginalChain)
    else:
        Upperbound = MaxAcceptantList[PositionOfLastOpenInPartialChain]
    return Upperbound


def CheckIfThisOneIsUnpairWillTheRemaingingFeasible(PartialChain, RNAOriginalChain, PositionOfLastOpenInPartialChain, PositionOfThisNucleotide, BPpiece):
    """
    ***NEW***
    In this function we are trying to find one feasible solution, under the assumption that 
    """
    RNAOriginalChain_copy = deepcopy(RNAOriginalChain)
    BPpiece_copy = deepcopy(BPpiece)
    RNAOriginalChain_copy[PositionOfThisNucleotide] = 'Unpaired'
    #print('The BPpiece Check:', BPpiece_copy)
    checkDic = {}
    LastTakenPosition = 0
    for i in range(len(BPpiece_copy)-1,-1,-1):
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
            if 'Taken'in RNAOriginalChain_copy:
                LastTakenPosition = len(RNAOriginalChain_copy) - RNAOriginalChain_copy[::-1].index('Taken') -1
            else:
                LastTakenPosition = 0
            #print('LastTakenPosition',LastTakenPosition)
            for j in range(max(PositionOfLastOpenInPartialChain+1+4,LastTakenPosition+1,PositionOfThisNucleotide+1), len(RNAOriginalChain_copy)):
                for Acceptant in AcceptantType:
                    if RNAOriginalChain[j] == Acceptant:
                        checkDic[str(i)] = 'paired'
                        RNAOriginalChain_copy[i] = 'TakenOpen'
                        RNAOriginalChain_copy[j] = 'Taken'
                        break
                else:
                    continue
                break
            #print('RNAchain now',RNAOriginalChain_copy)
        else:
            pass

    if checkDic == {}:
        check = True
    else:
        check = all(value == 'paired' for value in checkDic.values())
    #print('checkDic:',checkDic)
    #print('RNAchain now:',RNAOriginalChain_copy) 
    #print('check is:',check)		
    return check


def CheckIfThisPositionPreserved(PositionOfThisNucleotide, piece, BPpiece, PreservedList):
    """
    In this function, we check if in this position, the nucleotide is already preserved for pairing.
    NOTICE: we are not using this function now.
    """
    for item in PreservedList:
        if PositionOfThisNucleotide == item[2]:
            temp_dec = 1
            BPchoice = item[3]
            break
        else:
            temp_dec = 0
    if temp_dec == 1:
        result = BPchoice             
    else:
        result = 'Non-preserved'
    return result


def UpdatePotentialAcceptantListForOriginalChain():
    """
    In this function, we check the corresponding potential pair for the newly attached nucleotide.
    Also, when there is a i+k position has the same potential pair j as i, the list for i should be updated,
        the if j is the smallest of i's list, it will be removed from i's list to i+a's.
    NOTICE: we are not using this function now.
    """
    
    return

def UpdateBPChain(BPChain, ActionChosen):
    """
    This function updates the BP chain state after we attach a new one, it transfers uncomplete to newly complete: 1 to 1.5, 2 to 2.5
    """
    #print('test_ori_BP',BPChain)
    #NewBPChain = deepcopy(BPChain)
    NewBPChain = BPChain
    if ActionChosen.BasePair == 0:
        pass
    elif ActionChosen.BasePair == 1:
        #PositionOfAttachedNucleotide=ActionChosen.PositionOfNucleotide
        #PositionOfFirstCloseInPartialChain = NewBPChain.index(2)
        #NewBPChain[PositionOfAttachedNucleotide]=1.5
        #NewBPChain[PositionOfFirstCloseInPartialChain]=2.5
	pass
    elif ActionChosen.BasePair == 2:
        #print('FLAG ActionUpdate as 2')
        #print('Length of BPChain:',len(NewBPChain))
        #print('testNewBPChainOriginal:',NewBPChain)
        PositionOfAttachedNucleotide=ActionChosen.PositionOfNucleotide
        PositionOfLastOpenInPartialChain = len(NewBPChain) - 1 - NewBPChain[::-1].index(1)
        print('PositionOfAttachedNucleotide',PositionOfAttachedNucleotide)
        print('PositionOfLastOpenInPartialChain',PositionOfLastOpenInPartialChain)
        print('length of NewBPChain', len(NewBPChain))
        NewBPChain[PositionOfAttachedNucleotide]=2.5
        NewBPChain[PositionOfLastOpenInPartialChain]=1.5
    #print('testNewBPChainUpdated:',NewBPChain)
    return NewBPChain
    

def TakeAction(PartialChain, ActionChosen):
    """
    Here we take the action chosen and update the partial chain.
    $$$We also need to update BPchain.
    There are several situations:
        (1)Attach to left of the partial chain
        (2)Attach to right of the partial chain
    """
    #print TakeAction
    NewPartialChain = deepcopy(PartialChain)
    #NewPartialChain = PartialChain
    #print NewPartialChain.ChainItself
    #print NewPartialChain.BasePairChain
    NewPartialChain.ChainItself.append(ActionChosen.Nucleotide)
    NewPartialChain.BasePairChain.append(ActionChosen.BasePair)
    #print NewPartialChain.ChainItself
    #print NewPartialChain.BasePairChain
    #print('append',NewPartialChain.BasePairChain)
    NewPartialChain.BasePairChain = UpdateBPChain(NewPartialChain.BasePairChain,ActionChosen)
    #elif ActionChosen.AttachToDirection =='right':
        #NewPartialChain.ChainItself.insert(0, ActionChosen.Nucleotide)
        #NewPartialChain.BasePairChain.insert(0, ActionChosen.BasePair)
        ##print('append',NewPartialChain.BasePairChain)
        #NewPartialChain.BasePairChain = UpdateBPChain(PartialChain.BasePairChain,ActionChosen)
	#pass
    #print('NewPartialChain_Nucleotide', NewPartialChain.ChainItself, 'NewPartialChain_BasePairChain',NewPartialChain.BasePairChain)
    return NewPartialChain


def TransferIntoDP(PartialChain):
    """
    In this function, we transfer the BasePair notation into dot-and-bracket notation.
    """
    BPChain = PartialChain.BasePairChain
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



#$$$$$-----------------initialization--------------------:
"""
InitialPartialChain = PartialChain(ChainItself = [], ConnectivityLeft = 1, ConnectivityRight =1, BasePairChain = [], PreservedList = [], MaxAcceptantList = [None for x in range(len(RNAOriginalChain))])
PartialChainHistory = []
ActionHistory = []


print(InitialPartialChain.ChainItself)


while len(InitialPartialChain.ChainItself)<len(RNAOriginalChain):
    PartialChainHistory.append(InitialPartialChain.ChainItself)
    PossibleActionSet = GetPossibleActions(InitialPartialChain)
    print('PossibleActionSet is:',PossibleActionSet)
    ActionChosen = SelectFromPossibleActions(PossibleActionSet)
    ActionHistory.append(ActionChosen)
    print('ActionChosen:','Nucleotide:', ActionChosen.Nucleotide, 'AttachToDirection:', ActionChosen.AttachToDirection, 'AllowableOtherSide:', ActionChosen.AllowableOtherSide, 'BasePair:', ActionChosen.BasePair, 'PositionOfNucleotide:', ActionChosen.PositionOfNucleotide)
    InitialPartialChain = TakeAction(InitialPartialChain, ActionChosen)
    print('InitialPartialChain:', InitialPartialChain.ChainItself)
    print('InitialPartialChainBP:', InitialPartialChain.BasePairChain)
    
for item in ActionHistory:
    print('ActionChosen:','Nucleotide:', item.Nucleotide, 'AttachToDirection:', item.AttachToDirection, 'AllowableOtherSide:', item.AllowableOtherSide, 'BasePair:', item.BasePair, 'PositionOfNucleotide:', item.PositionOfNucleotide)

print(InitialPartialChain.BasePairChain)
print(TransferIntoDP(InitialPartialChain))

"""
