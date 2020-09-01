import numpy as np
import pandas as pd
import RNA

##------------------------------------------------------------------------------------------------------------------------------##
def IsPseudoFree(bp):
    bpstart = []
    nn = len(bp)
    for i in range(nn-1):
        if i+1 < bp[i]:
            bpstart.append(i+1)
    flag = 0
    bp_n = len(bpstart)
    for i in range(bp_n):
        for j in range(bp_n-i-1):
            if (bpstart[i]<bpstart[i+1+j] and bp[bpstart[i]-1] > bp[bpstart[i+1+j]-1]) or ( bp[bpstart[i]-1] < bpstart[i+1+j]):
                continue
            else:
                flag = 1
                break
    return flag
##################################################################################################################################


def HasBranch(DP):
    left_start = 0
    right_end = len(DP)
    for i in range(len(DP)):
        if DP[i] == "(":
            if i > left_start:
                left_start = i
        elif DP[i] == ")":
            if i < right_end:
                right_end = i


    if left_start < right_end:
        return 0
    else:
        return 1



##------------------------------------------------------------------------------------------------------------------------------##
def UnpairSegNum(segAssign):
####################################################
# 0: keep unpair
# 1: keep bp
# 2: all deleted
####################################################
    length = len(segAssign)
    segUnpair = np.zeros(length,dtype = int)
    Main_Start = np.amin(np.where(segAssign == 1))
    Main_End = np.amax(np.where(segAssign == 1))
    Sub_Start = np.amin(np.where(segAssign == 2))
    Sub_End = np.amax(np.where(segAssign == 2))
    for i in range(length):
        if segAssign[i] == 2:
            segUnpair[i] = 2
        elif segAssign[i] == 1:
            segUnpair[i] = 1
        elif segAssign[i] == 0:
            if i < Main_End and i > Main_Start:
                segUnpair[i] = 0
            elif i > Main_End and Main_End > Sub_End:
                segUnpair[i] = 2
            elif i < Main_Start and Main_Start < Sub_Start:
                segUnpair[i] = 2
            else:
                segUnpair[i] = 2
    return segUnpair
##################################################################################################################################




##------------------------------------------------------------------------------------------------------------------------------##
def ExtractSEQ(segBP,Seq):
    segSEQ = ""
    for i in range(len(segBP)):
        if segBP[i] >= 0:
            segSEQ += Seq[i]
    return segSEQ
##################################################################################################################################


##------------------------------------------------------------------------------------------------------------------------------##
def DeleteLeft(ii,segBP,arm,seg_i):
    rna_length = len(segBP)
    for i in range(rna_length):
        if i > ii and arm[i] == seg_i:
            segBP[i] = segBP[i] - 1
    return segBP
##################################################################################################################################


##------------------------------------------------------------------------------------------------------------------------------##
def DeleteMiddle(ii,segBP,fArray,arm,seg_i):
    rna_length = len(segBP)
    for i in range(rna_length):
        if fArray[i] - 1 > ii and arm[i] == seg_i:
            segBP[i] = segBP[i] - 1
    return segBP
##################################################################################################################################


##------------------------------------------------------------------------------------------------------------------------------##
def DeleteRight(ii,segBP,arm,seg_i):
    return segBP
##################################################################################################################################

def InsertAAAAA(segSEQFinal_Before,segDPFinal_Before):
    Insert_Start = 0
    Insert_End = len(segDPFinal_Before)
    for i in range(len(segDPFinal_Before)):
        if segDPFinal_Before[i] == "(":
            if i > Insert_Start:
                Insert_Start = i
        elif segDPFinal_Before[i] == ")":
            if i < Insert_End:
                Insert_End = i
    segSEQFinal =  segSEQFinal_Before[:Insert_Start+1]+"AAAAA"+segSEQFinal_Before[Insert_End:]
    segDPFinal = segDPFinal_Before[:Insert_Start+1]+"....."+segDPFinal_Before[Insert_End:]
    HairpinSeq = segSEQFinal_Before[Insert_Start]+"AAAAA"+segSEQFinal_Before[Insert_End]
    return HairpinSeq,segSEQFinal, segDPFinal






##------------------------------------------------------------------------------------------------------------------------------##
def  BP_to_DP(segBPFinal):
    segDP = ""
    for i in range(len(segBPFinal)):
        if segBPFinal[i] == 0:
            segDP += "."
        elif i + 1 < segBPFinal[i]:
            segDP += "("
        elif i + 1 > segBPFinal[i]:
            segDP += ")"
    return segDP
##########################################################################################################################################################################################
##########################################################################################################################################################################################
##########################################################################################################################################################################################


def Peudo_Decom(bp,seq,name):
    Name_Original = name
    BP_Original = bp
    Seq_Original = seq

    fArray = np.asarray(BP_Original,dtype = int)  # bp numpy array
    arm = np.ones(len(fArray),dtype=int) * -1
    armnum = 1
#####################################################################
    for nn in range(len(fArray)):
        if arm[nn] >= 0:
            continue
        elif fArray[nn] == 0:
            arm[nn] = 0
        elif nn == 0:
            arm[nn] = armnum
            arm[fArray[nn] - 1] = armnum
        else:
            if fArray[nn] == fArray[nn - 1] - 1:
                arm[nn] = armnum
                arm[fArray[nn] - 1] = armnum
            else:
                armnum += 1
                arm[nn] = armnum
                arm[fArray[nn] - 1] = armnum
##################### finish small segment assignment##################


#####################################################################
    unique, counts = np.unique(arm, return_counts=True)
    freq = dict(zip(unique, counts))
    freqList = sorted(freq, key = freq.get, reverse=True)
############# arm frequency list after sorting########################



########### combine arms into segments ###############
    segnum = np.ones(len(fArray),dtype=int) * -1
    segnum = 1
    segList = []
    segArray = np.zeros(len(fArray),dtype = int)
############### seg global info initialization ########################


#####################  Primary Chain############################
    segList.append(segnum)
    bpknotfreeMainArray = np.zeros(len(fArray), dtype=int)
    for i_1 in range(len(freqList)):
        if freqList[i_1] == 0:
            continue
        else:
            if IsPseudoFree(bpknotfreeMainArray + fArray * (arm == freqList[i_1])) == 0: # 0 means no pseudo
                bpknotfreeMainArray += fArray * (arm == freqList[i_1])
                segArray += segnum * (arm == freqList[i_1])
                freqList[i_1] = 0
###############################################################



    while sum(freqList) > 0:
        bpknotfreeArray = np.zeros(len(fArray), dtype=int)
        segnum += 1
        segList.append(segnum)
        # segStartIdx = 0
        # segEndIdx = 0
        for i_22 in range(len(freqList)):
            if freqList[i_22] > 0:
                bpknotfreeArray += fArray * (arm == freqList[i_22])
                segArray += segnum * (arm == freqList[i_22])
                freqList[i_22] = 0
                break

        for i_2 in range(len(freqList)):
            if freqList[i_2] > 0 and IsPseudoFree(bpknotfreeArray + fArray * (arm == freqList[i_2])) == 0:
                TempDP = BP_to_DP(bpknotfreeArray + fArray * (arm == freqList[i_2]))
                left_start = TempDP.index("(")
                left_end = len(fArray)
                right_start = TempDP.index(")")
                right_end = 0
                for ii in range(len(fArray)):
                    if TempDP[ii] == "(":
                        left_end = ii
                    elif TempDP[ii] == ")":
                        right_end = ii
                if HasBranch(TempDP) == 0 and sum(bpknotfreeMainArray[left_start:left_end+1]) == 0 and sum(bpknotfreeMainArray[right_start:right_end+1]) == 0:
                    bpknotfreeArray += fArray * (arm == freqList[i_2])
                    segArray += segnum * (arm == freqList[i_2])
                    freqList[i_2] = 0

    kfe = 0
    bfe = 0
    for i_3 in segList:
        if i_3 == 1:
            segBPFinal = fArray*(segArray == i_3)
            segSEQFinal = Seq_Original.replace("\n","")
            segDPFinal = BP_to_DP(segBPFinal)

            b_temp = RNA.fold_compound(segSEQFinal)
            b_temp.pf()
            bfe = b_temp.eval_structure(segDPFinal)

        else:
            
            ####################################################
            # 0: Unpaired
            # 1: Main
            # 2: Paired but not in this structure
            ####################################################
            segNew = 2*(segArray != i_3)-2*(segArray == 0)+1*(segArray == i_3)
            segUnpair = UnpairSegNum(segNew)
            segArm = np.delete(np.unique(arm*(segArray == i_3)),0)
            segBP = np.copy(fArray)
            rna_length = len(fArray)
            for seg_i in segArm:
                seg_Start = np.amin(np.where(arm == seg_i))
                seg_End = np.amax(np.where(arm == seg_i))
                for ii in range(rna_length):
                   if segUnpair[ii] == 2:
                        segBP[ii] = -10
                        if ii < seg_Start:
                            segBP = DeleteLeft(ii,segBP,arm,seg_i)
                        elif seg_Start < ii and ii < seg_End:
                            segBP = DeleteMiddle(ii,segBP,fArray,arm,seg_i)
                        elif ii > seg_End:
                            segBP = DeleteRight(ii,segBP,arm,seg_i)
                   else:
                       continue
            idxs = (segBP != -10)
            segBPFinal = segBP[idxs]
            segSEQFinal_Before = ExtractSEQ(segBP,Seq_Original.replace("\n",""))
            segDPFinal_Before = BP_to_DP(segBPFinal)
            HairpinSeq,segSEQFinal,segDPFinal = InsertAAAAA(segSEQFinal_Before,segDPFinal_Before)
            
            k_temp = RNA.fold_compound(segSEQFinal)
            k_temp.pf()
            kfe += k_temp.eval_structure(segDPFinal)
    return bfe, kfe



