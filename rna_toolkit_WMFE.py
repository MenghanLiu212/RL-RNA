import numpy as np
import math



def is_pseudoknotted(bp):
    bp_paired_list = []
    for i in range(len(bp)):
        if i < bp[i]:
            bp_paired_list.append([i,bp[i]-1])
    is_pseudoknotted_idx = 0
    for ii in range(len(bp_paired_list)):
        test_paired = bp_paired_list[ii]
        for jj in range(len(bp_paired_list)):
            if test_paired[0] < bp_paired_list[jj][0] and test_paired[1] > bp_paired_list[jj][0] and test_paired[1] < bp_paired_list[jj][1]:
                is_pseudoknotted_idx += 1
    return is_pseudoknotted_idx

def bp_to_dp(bp):
    dp = ''
    for i in range(len(bp)):
        if bp[i] == 0:
            dp += '.'
        elif i < bp[i]:
            dp += '('
        else:
            dp += ')'
    return dp


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



def entropy(seq, mm):
    seqDict = {}
    n = len(seq)
    for i in range(n-mm+1):
        if seq[i:i+mm] in seqDict:
            seqDict[seq[i:i+mm]] += 1
        else:
            seqDict[seq[i:i+mm]] = 1
    nn = float(n-mm+1)
    x = 0
    for i in seqDict:
        prob = seqDict[i]/nn
        seqDict[i] = - prob * math.log(prob, 2)
        x += seqDict[i]
    return x


def entropy_max(nn, ent_n):
    entropy_max_value = 0
    if nn <= ent_n:
        entropy_max_value = 0
    elif nn - ent_n + 1 <= 4**ent_n:
            prob = 1/float(nn - ent_n + 1)
            for i in range(nn - ent_n + 1):
                entropy_max_value +=  (-(prob) * math.log(prob, 2))
    elif nn - ent_n + 1 > 4**ent_n:
        a = (nn - ent_n + 1) / (4**ent_n)
        b = (nn - ent_n + 1) % (4**ent_n)
        entropy_max_value = (-b) * (a+1)/float(nn - ent_n + 1) * math.log((a+1)/float(nn - ent_n + 1), 2) - ((4**ent_n-b) * a / float(nn - ent_n + 1) * math.log(a/float(nn - ent_n + 1), 2))
    return entropy_max_value