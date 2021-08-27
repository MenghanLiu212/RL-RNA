import numpy as np
from re import finditer
from os import listdir
import matplotlib.pyplot as plt
import matplotlib as mpl
from sklearn.linear_model import LinearRegression
from scipy.optimize import curve_fit
from scipy.stats import ttest_ind

def calc_MCC(prediction, actual, name, out_number, alg):
    TP = 0
    TN = 0
    FP = 0
    FN = 0
    with open("{}_{}_{}_colormap.col".format(name, out_number, alg), 'w+') as f:

        for i, (a, b) in enumerate(zip(prediction, actual)):
            #true
            if a == b:
                f.write("{}:green ".format(i + 1)) #forna is 1-indexed...
                #negative
                if b == -1:
                    TN += 1
                #positive
                else:
                    TP += 1

            #false
            else:
                f.write("{}:red ".format(i + 1))
                #positive
                if b == -1:
                    FP += 1
                #negative
                else:
                    FN += 1

    MCC = ((TP * TN) - (FP * FN)) / np.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

    return MCC

#I represent structures as a pairing vector
def parse_dot_bracket(input):
    output = np.full(len(input), -1)
    #I'm not sure this is the most efficent way to do this, but I'm lazy.
    more = True
    while more:
        more = False

        #finds matched parenthesis
        for x in finditer(r"\([^()]*\)", input):
            more = True
            output[x.start()] = x.end()-1
            output[x.end()-1] = x.start()

            #its recursive...
            input=input[0:x.start()] + "." + input[x.start()+1:x.end()-1] + "." + input[x.end():]

    return output

def make_fasta_db(name, out_number, alg, seq, struct):
    with open("{}_{}_{}.fasta".format(name, out_number, alg), 'w+') as f:
        f.write(">"+name+"\n")
        f.write(seq+"\n")
        f.write(struct+"\n")

def main(rank):

    step = 4
    out_number = rank
    if rank == 0:
        rank = 1
        step = 1

    #read all files in the directory, extract dot-brackets and calculate MCC
    data_fold = []
    data_roll = []
    data_diff = []
    actual_foldabilityMFE = []
    RNAfold_foldabilityMFE = []
    rollout_foldabilityMFE = []
    actual_foldabilityNFE = []
    RNAfold_foldabilityNFE = []
    rollout_foldabilityNFE = []
    runtime = []
    actual_MFE = []
    RNAfold_MFE = []
    rollout_MFE = []
    names = []
    seq_lens = []
    families = []
    count = 0
    for fname in listdir():
        if fname.split(".")[-1] == "csv":
            with open(fname, "r") as f:
                lines = f.readlines()
                for l in lines[rank::step]: #first line is a header, so ranks are 1-indexed.
                    l = l.split(",")
                    if l[0] == 'NONE':
                        continue
                    name = l[0]
                    try:
                        seq = l[1]
                    except:
                        print(l)
                    seq_lens.append(len(seq))
                    actual = l[2]
                    RNAfold = l[3]
                    try:
                        rollout = l[4]
                    except:
                        print(l)
                    family = name.split("_")[0]
                    actual_foldabilityMFE.append(float(l[13]))
                    RNAfold_foldabilityMFE.append(float(l[15]))
                    rollout_foldabilityMFE.append(float(l[17]))
                    actual_foldabilityNFE.append(float(l[14]))
                    RNAfold_foldabilityNFE.append(float(l[16]))
                    rollout_foldabilityNFE.append(float(l[18]))
                    if rank == 1 and step == 4 and float(l[16]) > float(l[18]) and 'NFE' in fname and not 'MFE' in fname:
                        count += 1
                        print(count, name, float(l[16]), float(l[18])) #structures where RNAfold had higher foldability than ExpertRNA
                    elif rank == 1 and step ==4 and float(l[15]) > float(l[17]) and 'MFE' in fname and not 'NFE' in fname:
                        count += 1
                        print(count, name, float(l[15]), float(l[17])) #structures where RNAfold had higher foldability than expertRNA
                    runtime.append(float(l[12]))
                    actual_MFE.append(float(l[19]))
                    RNAfold_MFE.append(float(l[20]))
                    rollout_MFE.append(float(l[21]))
                    
                    make_fasta_db(name, out_number, "actual", seq, actual)
                    make_fasta_db(name, out_number, "RNAfold", seq, RNAfold)
                    make_fasta_db(name, out_number, "rollout", seq, rollout)
                    actual = parse_dot_bracket(actual)            
                    rollout = parse_dot_bracket(rollout)
                    RNAfold = parse_dot_bracket(RNAfold)
                    acc_fold = calc_MCC(RNAfold, actual, name, out_number, "RNAfold")
                    acc_roll = calc_MCC(rollout, actual, name, out_number, "rollout")
                    data_fold.append(acc_fold)
                    data_roll.append(acc_roll)
                    data_diff.append(acc_roll - acc_fold)
                    names.append(name)
                    families.append(family)
                    #print(name, acc_fold, acc_roll, acc_roll-acc_fold)

    all_data = np.array([data_fold, data_roll, data_diff, actual_foldabilityMFE, RNAfold_foldabilityMFE, rollout_foldabilityMFE, seq_lens, runtime, actual_MFE, RNAfold_MFE, rollout_MFE, actual_foldabilityNFE, RNAfold_foldabilityNFE, rollout_foldabilityNFE], dtype=float)
    DATA_FOLD = 0
    DATA_ROLL = 1
    DATA_DIFF = 2
    ACTUAL_FOLDABILITY_MFE = 3    
    RNAFOLD_FOLDABILITY_MFE = 4
    ROLLOUT_FOLDABILITY_MFE = 5
    SEQ_LEN = 6
    RUNTIME = 7
    ACTUAL_MFE = 8
    RNAFOLD_MFE = 9
    ROLLOUT_MFE = 10
    ACTUAL_FOLDABILITY_NFE = 11
    RNAFOLD_FOLDABILITY_NFE = 12
    ROLLOUT_FOLDABILITY_NFE = 13
    sort = all_data[DATA_DIFF].argsort()
    all_data = all_data[:,sort]
    names = np.array(names, dtype=str)[sort]
    families = np.array(families)[sort]
    good_MFE = all_data[ACTUAL_MFE] < 5

    statistic, pvalue = ttest_ind(all_data[DATA_FOLD], all_data[DATA_ROLL])

    print("{},{:.4f},{:.4f},{:.4f},{:.4f},{:.4f}".format(out_number, np.mean(all_data[DATA_ROLL]), np.median(all_data[DATA_ROLL]), np.mean(all_data[DATA_DIFF]), np.std(all_data[DATA_DIFF]), pvalue))

    with open("scores_{}.dat".format(out_number), 'w') as f:
        for n, d in zip(names, all_data.T):
            f.write("{} {} {} {}\n".format(n, d[DATA_FOLD], d[DATA_ROLL], d[DATA_DIFF]))


    mpl.rcParams.update({'font.size': 12, 'lines.markersize' : 2})
    mpl.rc('text', usetex=True)
    mpl.rc('xtick', labelsize=8)
    mpl.rc('ytick', labelsize=8)
    W2 = np.array(plt.rcParams["figure.figsize"])*2.6/plt.rcParams["figure.figsize"][0] #for 2 graphs across with the same width as the text column in PLOS
    W3 = np.array(plt.rcParams["figure.figsize"])*2.5/plt.rcParams["figure.figsize"][0] #for 3 graphs across the whole page in PLOS
    WeqH = np.array([2.6, 2.6])

    fnames = sorted(list(set(families)))
    if len(fnames) == 5: #rfam dataset has 5 families, matthews has 3, make the colors line up
        fnames[1], fnames[2], fnames[3] = fnames[2], fnames[3], fnames[1]
    best_RNAfold = all_data[DATA_FOLD] > 0.9
    best_rollout = all_data[DATA_ROLL] > 0.9
    """   
    less_best_rollout = all_data[DATA_ROLL] > 0.8
    less_worst_RNAfold = all_data[DATA_FOLD] < 0.8
    sameMCC = all_data[DATA_ROLL] == all_data[DATA_FOLD]
    good_roll_bad_fold = names[~sameMCC & less_best_rollout & less_worst_RNAfold]
    print("number not same MCC", len(names[~sameMCC]))
    print("number rolloutMCC > 0.8 and RNAfoldMCC < 0.8: ", len(good_roll_bad_fold))
    print("as a percentage: ", len(good_roll_bad_fold) / len(names[~sameMCC]))
    FEdifference = (all_data[ROLLOUT_MFE] - all_data[RNAFOLD_MFE]) / all_data[RNAFOLD_MFE]
    FEdifference5 = FEdifference > 0.05
    good_roll_bad_fold_big_diff = names[~sameMCC & less_best_rollout & less_worst_RNAfold & FEdifference5]
    print("number with FE diff (roll-fold/fold) > 5%: ", len(good_roll_bad_fold_big_diff))
    print("as a percentage: ", len(good_roll_bad_fold_big_diff) / len(names[~sameMCC]))
    foldability_difference = (all_data[ROLLOUT_FOLDABILITY_NFE] - all_data[RNAFOLD_FOLDABILITY_NFE]) / all_data[RNAFOLD_FOLDABILITY_NFE]
    foldability_difference5 = foldability_difference > 0.05
    good_roll_bad_fold_big_diff = names[~sameMCC & less_best_rollout & less_worst_RNAfold & foldability_difference5]
    print("number with foldability diff > 5%: ", len(good_roll_bad_fold_big_diff))
    print("as a percentage: ", len(good_roll_bad_fold_big_diff) / len(names[~sameMCC]))

    print()


    less_worst_rollout = all_data[DATA_ROLL] < 0.8
    less_best_RNAfold = all_data[DATA_FOLD] > 0.8
    sameMCC = all_data[DATA_ROLL] == all_data[DATA_FOLD]
    good_fold_bad_roll = names[~sameMCC & less_worst_rollout & less_best_RNAfold]
    print("number RNAfoldMCC > 0.8 and rolloutMCC < 0.8: ", len(good_fold_bad_roll))
    print("as a percentage: ", len(good_fold_bad_roll) / len(names[~sameMCC]))
    FEdifference = (all_data[ROLLOUT_MFE] - all_data[RNAFOLD_MFE]) / all_data[RNAFOLD_MFE]
    FEdifference5 = FEdifference < -0.05
    good_fold_bad_roll_big_diff = names[~sameMCC & less_worst_rollout & less_best_RNAfold & FEdifference5]
    print("number with FE diff (roll-fold/fold) < -5%: ", len(good_fold_bad_roll_big_diff))
    print("as a percentage: ", len(good_fold_bad_roll_big_diff) / len(names[~sameMCC]))
    foldability_difference = (all_data[ROLLOUT_FOLDABILITY_NFE] - all_data[RNAFOLD_FOLDABILITY_NFE]) / all_data[RNAFOLD_FOLDABILITY_NFE]
    foldability_difference5 = foldability_difference > 0.05
    good_fold_bad_roll_big_diff = names[~sameMCC & less_worst_rollout & less_best_RNAfold & foldability_difference5]
    print("number with foldability diff > 5%: ", len(good_fold_bad_roll_big_diff))
    print("as a percentage: ", len(good_fold_bad_roll_big_diff) / len(names[~sameMCC]))
    """
    #comparison scatterplot
    fig, ax = plt.subplots(figsize=WeqH)
    for fname in fnames:
        ax.scatter(all_data[DATA_FOLD][families == fname], all_data[DATA_ROLL][families == fname], label=fname, alpha=0.4)
    line = np.linspace(-1, 1, 20)
    ax.plot(line, line, c='k', linewidth=0.75)
    ax.set_xlabel("RNAfold MCC")
    ax.set_ylabel("ExpertRNA MCC")
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.1, 1.1)
    #ax.legend(fontsize=8), bbox_to_anchor=(1.04, 1))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)  
    ax.axis('scaled')  
    plt.tight_layout()
    plt.savefig("scatter_comparison_{}.png".format(out_number), dpi=600)
    plt.close()

    #histograms!
    #overlay
    fig, ax = plt.subplots(figsize=W2)
    ax.hist(data_fold, bins=50, range=(-1, 1), alpha=0.4, label="RNAfold")
    ax.hist(data_roll, bins=50, range=(-1, 1), alpha=0.4, label="ExpertRNA")
    ax.set_xlabel("MCC")
    ax.set_ylabel("Number")
    ax.set_xlim(-1, 1)
    ax.legend(fontsize=8)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.tight_layout()
    plt.savefig("comparison_{}.png".format(out_number), dpi=600)
    plt.close()

    #difference
    fig, ax = plt.subplots(figsize=W2)
    ax.hist(data_diff, range=(-1.4,1.4), bins=28)
    ax.set_xlabel("MCC difference")
    ax.set_ylabel("Number")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlim(-1.4, 1.4)
    plt.tight_layout()
    plt.savefig("diff_{}.png".format(out_number), dpi=600)
    plt.close()

    #runtime
    fig, ax = plt.subplots(figsize=W2)
    def func(x, a, b, c): #exponential fit
        return a * np.power(2, b*x) + c
    def func2(x, a, b, c): #polynomial fit (worked better)
        return a*np.power(x, 2) + b * x + c
    #popt, pcov = curve_fit(func, all_data[SEQ_LEN], all_data[RUNTIME])
    popt2, pcov2 = curve_fit(func2, all_data[SEQ_LEN], all_data[RUNTIME])
    mi = min(all_data[SEQ_LEN])
    ma = max(all_data[SEQ_LEN])
    modelX = np.linspace(mi, ma, 100)
    modelY = func2(modelX, *popt2)
    ax.scatter(all_data[SEQ_LEN], all_data[RUNTIME], alpha=0.4)
    ax.plot(modelX, modelY, label=r"quadratic fit ${:.2f}x^{{2}} + {:.2f}x + {:.2f}$".format(popt2[0], popt2[1], popt2[2]), linewidth=0.75, c='k')
    ax.set_xlabel("Sequence Length")
    ax.set_ylabel("Runtime (s)")
    ax.legend(fontsize=4)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.tight_layout()
    plt.savefig("runtime_{}.png".format(out_number), dpi=600)
    plt.close()

    max_mfe = max([max(all_data[ACTUAL_MFE][good_MFE] / all_data[SEQ_LEN][good_MFE]), max(all_data[ROLLOUT_MFE][good_MFE] / all_data[SEQ_LEN][good_MFE]), max(all_data[RNAFOLD_MFE][good_MFE] / all_data[SEQ_LEN][good_MFE])]) + 0.02
    min_mfe = min([min(all_data[ACTUAL_MFE][good_MFE] / all_data[SEQ_LEN][good_MFE]), min(all_data[ROLLOUT_MFE][good_MFE] / all_data[SEQ_LEN][good_MFE]), min(all_data[RNAFOLD_MFE][good_MFE] / all_data[SEQ_LEN][good_MFE])]) - 0.02

    #actual vs rollout MFE
    fig, ax = plt.subplots(figsize=WeqH)
    line = np.linspace(-0.7, 0, 20)
    for fname in fnames:
        ax.scatter(all_data[ACTUAL_MFE][good_MFE][families[good_MFE] == fname] / all_data[SEQ_LEN][good_MFE][families[good_MFE] == fname], all_data[ROLLOUT_MFE][good_MFE][families[good_MFE] == fname] / all_data[SEQ_LEN][good_MFE][families[good_MFE] == fname], label=fname, alpha=0.4)
    ax.scatter(all_data[ACTUAL_MFE][best_rollout] / all_data[SEQ_LEN][best_rollout], all_data[ROLLOUT_MFE][best_rollout] / all_data[SEQ_LEN][best_rollout], c='red', s=1, label="Best ExpertRNA")
    ax.scatter(all_data[ACTUAL_MFE][best_RNAfold] / all_data[SEQ_LEN][best_RNAfold], all_data[ROLLOUT_MFE][best_RNAfold] / all_data[SEQ_LEN][best_RNAfold], c='cyan', s=0.5, label="Best RNAfold")
    ax.plot(line, line, c='k', linewidth=0.75)
    ax.set_xlabel("Actual FE")
    ax.set_ylabel("ExpertRNA FE")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlim(min_mfe, max_mfe)
    ax.set_ylim(min_mfe, max_mfe)
    #ax.axis('scaled')
    #leg = ax.legend(fontsize=8, bbox_to_anchor=(1.04, 1))
    plt.tight_layout()
    plt.savefig("actual_roll_MFE_comparison_{}.png".format(out_number), dpi=600)#, bbox_extra_artists=[leg], bbox_inches='tight')
    plt.close()

    #actual vs RNAfold MFE
    fig, ax = plt.subplots(figsize=WeqH)
    for fname in fnames:
        ax.scatter(all_data[ACTUAL_MFE][good_MFE][families[good_MFE] == fname] / all_data[SEQ_LEN][good_MFE][families[good_MFE] == fname], all_data[RNAFOLD_MFE][good_MFE][families[good_MFE] == fname] / all_data[SEQ_LEN][good_MFE][families[good_MFE] == fname], label=fname, alpha=0.4)
    ax.scatter(all_data[ACTUAL_MFE][best_rollout] / all_data[SEQ_LEN][best_rollout], all_data[RNAFOLD_MFE][best_rollout] / all_data[SEQ_LEN][best_rollout], c='red', s=1, label="Best ExpertRNA")
    ax.scatter(all_data[ACTUAL_MFE][best_RNAfold] / all_data[SEQ_LEN][best_RNAfold], all_data[RNAFOLD_MFE][best_RNAfold] / all_data[SEQ_LEN][best_RNAfold], c='cyan', s=0.5, label="Best RNAfold")
    ax.plot(line, line, c='k', linewidth=0.75)
    ax.set_xlabel("Actual FE")
    ax.set_ylabel("RNAfold FE")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlim(min_mfe, max_mfe)
    ax.set_ylim(min_mfe, max_mfe)
    #ax.axis('scaled')
    #leg = ax.legend(fontsize=8, bbox_to_anchor=(1.04, 1))
    plt.tight_layout()
    plt.savefig("actual_RNAfold_MFE_comparison_{}.png".format(out_number), dpi=600)#, bbox_extra_artists=[leg], bbox_inches='tight')
    plt.close()

    #RNAfold MFE vs Rollout MFE
    fig, ax = plt.subplots(figsize=WeqH)
    #line = np.linspace(-0.7, -0.1, 20)
    for fname in fnames:
        ax.scatter(all_data[RNAFOLD_MFE][families == fname] / all_data[SEQ_LEN][families == fname], all_data[ROLLOUT_MFE][families == fname] / all_data[SEQ_LEN][families == fname], label=fname, alpha=0.4)
    ax.scatter(all_data[RNAFOLD_MFE][best_rollout] / all_data[SEQ_LEN][best_rollout], all_data[ROLLOUT_MFE][best_rollout] / all_data[SEQ_LEN][best_rollout], c='red', s=1, label="Best ExpertRNA")
    ax.scatter(all_data[RNAFOLD_MFE][best_RNAfold] / all_data[SEQ_LEN][best_RNAfold], all_data[ROLLOUT_MFE][best_RNAfold] / all_data[SEQ_LEN][best_RNAfold], c='cyan', s=0.5, label="Best RNAfold")
    ax.plot(line, line, c='k', linewidth=0.75)
    ax.set_xlabel("RNAfold FE")
    ax.set_ylabel("ExpertRNA FE")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlim(min_mfe, max_mfe)
    ax.set_ylim(min_mfe, max_mfe)
    #ax.axis('scaled')
    #leg = ax.legend(fontsize=8, bbox_to_anchor=(1.04, 1))
    fig.tight_layout()
    plt.savefig("RNAfold_roll_MFE_comparison_{}.png".format(out_number), dpi=600)#, bbox_extra_artists=[leg], bbox_inches='tight')
    plt.close()

    #actual vs rollout foldability
    fig, ax = plt.subplots(figsize=WeqH)
    line = np.linspace(0, 1, 20)
    for fname in fnames:
        ax.scatter(all_data[ACTUAL_FOLDABILITY_MFE][families == fname], all_data[ROLLOUT_FOLDABILITY_MFE][families == fname], label=fname, alpha=0.4)
    ax.scatter(all_data[ACTUAL_FOLDABILITY_MFE][best_rollout], all_data[ROLLOUT_FOLDABILITY_MFE][best_rollout], c='red', s=1, label="Best ExpertRNA")
    ax.scatter(all_data[ACTUAL_FOLDABILITY_MFE][best_RNAfold], all_data[ROLLOUT_FOLDABILITY_MFE][best_RNAfold], c='cyan', s=0.5, label="Best RNAfold")
    ax.plot(line, line, c='k', linewidth=0.75)
    ax.set_xlabel("Actual foldability", wrap=True)
    ax.set_ylabel("ExpertRNA foldability", wrap=True)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.axis('scaled')
    #leg = ax.legend(fontsize=8, bbox_to_anchor=(1.04, 1))
    fig.tight_layout()
    plt.savefig("actual_roll_foldabilityMFE_comparison_{}.png".format(out_number), dpi=600)#, bbox_extra_artists=[leg], bbox_inches='tight')
    plt.close()

    #actual vs RNAfold foldability
    fig, ax = plt.subplots(figsize=WeqH)
    line = np.linspace(0, 1, 20)
    for fname in fnames:
        ax.scatter(all_data[ACTUAL_FOLDABILITY_MFE][families == fname], all_data[RNAFOLD_FOLDABILITY_MFE][families == fname], label=fname, alpha=0.4)
    ax.scatter(all_data[ACTUAL_FOLDABILITY_MFE][best_rollout], all_data[RNAFOLD_FOLDABILITY_MFE][best_rollout], c='red', s=1, label="Best ExpertRNA")
    ax.scatter(all_data[ACTUAL_FOLDABILITY_MFE][best_RNAfold], all_data[RNAFOLD_FOLDABILITY_MFE][best_RNAfold], c='cyan', s=0.5, label="Best RNAfold")
    ax.plot(line, line, c='k', linewidth=0.75)
    ax.set_xlabel("Actual foldability")
    ax.set_ylabel("RNAfold foldability")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.axis('scaled')
    #leg = ax.legend(fontsize=8, bbox_to_anchor=(1.04, 1))
    fig.tight_layout()
    plt.savefig("actual_RNAfold_foldabilityMFE_comparison_{}.png".format(out_number), dpi=600)#, bbox_extra_artists=[leg], bbox_inches='tight')
    plt.close()

    #RNAfold vs rollout foldability
    fig, ax = plt.subplots(figsize=WeqH)
    line = np.linspace(0, 1, 20)
    for fname in fnames:
        ax.scatter(all_data[RNAFOLD_FOLDABILITY_MFE][families == fname], all_data[ROLLOUT_FOLDABILITY_MFE][families == fname], label=fname, alpha=0.4)
    ax.scatter(all_data[RNAFOLD_FOLDABILITY_MFE][best_rollout], all_data[ROLLOUT_FOLDABILITY_MFE][best_rollout], c='red', s=1, label="Best ExpertRNA")
    ax.scatter(all_data[RNAFOLD_FOLDABILITY_MFE][best_RNAfold], all_data[ROLLOUT_FOLDABILITY_MFE][best_RNAfold], c='cyan', s=0.5, label="Best RNAfold")
    ax.plot(line, line, c='k', linewidth=0.75)
    ax.set_xlabel("RNAfold foldability")
    ax.set_ylabel("ExpertRNA foldability")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.axis('scaled')
    #leg = ax.legend(fontsize=8, bbox_to_anchor=(1.04, 1))
    fig.tight_layout()
    plt.savefig("RNAfold_roll_foldabilityMFE_comparison_{}.png".format(out_number), dpi=600)#, bbox_extra_artists=[leg], bbox_inches='tight')
    plt.close()

    #actual vs rollout foldability
    fig, ax = plt.subplots(figsize=WeqH)
    line = np.linspace(0, 1, 20)
    for fname in fnames:
        ax.scatter(all_data[ACTUAL_FOLDABILITY_NFE][families == fname], all_data[ROLLOUT_FOLDABILITY_NFE][families == fname], label=fname, alpha=0.4)
    ax.scatter(all_data[ACTUAL_FOLDABILITY_NFE][best_rollout], all_data[ROLLOUT_FOLDABILITY_NFE][best_rollout], c='red', s=1, label="Best ExpertRNA")
    ax.scatter(all_data[ACTUAL_FOLDABILITY_NFE][best_RNAfold], all_data[ROLLOUT_FOLDABILITY_NFE][best_RNAfold], c='cyan', s=0.5, label="Best RNAfold")
    ax.plot(line, line, c='k', linewidth=0.75)
    ax.set_xlabel("Actual foldability")
    ax.set_ylabel("ExpertRNA foldability")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.axis('scaled')
    #leg = ax.legend(fontsize=8, bbox_to_anchor=(1.04, 1))
    fig.tight_layout()
    plt.savefig("actual_roll_foldabilityNFE_comparison_{}.png".format(out_number), dpi=600)#, bbox_extra_artists=[leg], bbox_inches='tight')
    plt.close()

    #actual vs RNAfold foldability
    fig, ax = plt.subplots(figsize=WeqH)
    line = np.linspace(0, 1, 20)
    for fname in fnames:
        ax.scatter(all_data[ACTUAL_FOLDABILITY_NFE][families == fname], all_data[RNAFOLD_FOLDABILITY_NFE][families == fname], label=fname, alpha=0.4)
    ax.scatter(all_data[ACTUAL_FOLDABILITY_NFE][best_rollout], all_data[RNAFOLD_FOLDABILITY_NFE][best_rollout], c='red', s=1, label="Best ExpertRNA")
    ax.scatter(all_data[ACTUAL_FOLDABILITY_NFE][best_RNAfold], all_data[RNAFOLD_FOLDABILITY_NFE][best_RNAfold], c='cyan', s=0.5, label="Best RNAfold")
    ax.plot(line, line, c='k', linewidth=0.75)
    ax.set_xlabel("Actual foldability")
    ax.set_ylabel("RNAfold foldability")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.axis('scaled')
    #leg = ax.legend(fontsize=8, bbox_to_anchor=(1.04, 1))
    fig.tight_layout()
    plt.savefig("actual_RNAfold_foldabilityNFE_comparison_{}.png".format(out_number), dpi=600)#, bbox_extra_artists=[leg], bbox_inches='tight')
    plt.close()

    #RNAfold vs rollout foldability
    fig, ax = plt.subplots(figsize=WeqH)
    line = np.linspace(0, 1, 20)
    for fname in fnames:
        ax.scatter(all_data[RNAFOLD_FOLDABILITY_NFE][families == fname], all_data[ROLLOUT_FOLDABILITY_NFE][families == fname], label=fname, alpha=0.4)
    ax.scatter(all_data[RNAFOLD_FOLDABILITY_NFE][best_rollout], all_data[ROLLOUT_FOLDABILITY_NFE][best_rollout], c='red', s=1, label="Best ExpertRNA")
    ax.scatter(all_data[RNAFOLD_FOLDABILITY_NFE][best_RNAfold], all_data[ROLLOUT_FOLDABILITY_NFE][best_RNAfold], c='cyan', s=0.5, label="Best RNAfold")
    ax.plot(line, line, c='k', linewidth=0.75)
    ax.set_xlabel("RNAfold foldability")
    ax.set_ylabel("ExpertRNA foldability")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.axis('scaled')
    fig.tight_layout()
    plt.savefig("RNAfold_roll_foldabilityNFE_comparison_{}.png".format(out_number), dpi=600)#, bbox_extra_artists=[leg], bbox_inches='tight')
    plt.close()

    #separate legend
    leg = ax.legend(fontsize=12, framealpha=1)
    fig2 = leg.figure
    fig2.canvas.draw()
    bbox = leg.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig2.savefig("legend.png", dpi="figure", bbox_inches=bbox)

    #how accurate is foldability??
    #plt.figure()
    #for fname in fnames:
    #    model = LinearRegression()
    #    model.fit(all_data[ROLLOUT_FOLDABILITY][families == fname].reshape(-1, 1), all_data[DATA_ROLL][families == fname].reshape(-1, 1))
    #    mi = min(all_data[ROLLOUT_FOLDABILITY])
    #    ma = max(all_data[ROLLOUT_FOLDABILITY])
    #    modelX =  np.linspace(mi, ma, 100)
    #    modelY = model.predict(modelX[:, np.newaxis])
    #    plt.scatter(all_data[ROLLOUT_FOLDABILITY][families == fname], all_data[DATA_ROLL][families == fname], label=fname+' {:.2f}'.format(model.score(all_data[ROLLOUT_FOLDABILITY][families == fname].reshape(-1, 1), all_data[DATA_ROLL][families == fname].reshape(-1, 1))))
    #    plt.plot(modelX, modelY)
    #plt.xlabel("ENTRNA foldability")
    #plt.ylabel("MCC")
    #plt.legend(fontsize=14)
    #plt.tight_layout()
    #plt.savefig("foldability_{}.png".format(out_number))
    #plt.close()

    #plt.figure(figsize=(len(all_data[DATA_FOLD])/18, 20))
    #ma = np.max(np.array([np.max(all_data[RNAFOLD_FOLDABILITY]), np.max(all_data[ROLLOUT_FOLDABILITY])]))
    #mi = np.min(np.array([np.min(all_data[RNAFOLD_FOLDABILITY]), np.min(all_data[ROLLOUT_FOLDABILITY])]))
    ##cmap = mpl.cm.viridis
    ##norm = mpl.colors.Normalize(vmin=mi, vmax=ma)
    #a = plt.scatter(range(len(all_data[DATA_FOLD])), all_data[DATA_FOLD], c=all_data[RNAFOLD_FOLDABILITY], cmap='viridis', vmin=mi, vmax=ma, marker='o', label="RNAfold")
    #b = plt.scatter(range(len(all_data[DATA_FOLD])), all_data[DATA_ROLL], c=all_data[ROLLOUT_FOLDABILITY], cmap='viridis', vmin=mi, vmax=ma, marker='v', label="Rollout")
    #plt.legend(handles=[a, b], fontsize=14)
    #plt.colorbar(label='Free energy')
    #plt.xlabel("Structure")
    #plt.ylabel("MCC")
    #plt.xticks(range(len(all_data[DATA_FOLD])), [n.replace('_', r'\_') for n in names], rotation='vertical', fontsize=6)
    #plt.xlim((-0.5, len(all_data[DATA_FOLD])+0.5))
    #plt.ylim((-1.1, 1.1))
    #plt.tight_layout() 
    #plt.savefig("all_data_{}.png".format(out_number))
    #plt.close()
    

if __name__ == "__main__":
    print('branch,avg_mcc,median_mcc,avg_improvement,stdev_improvement,pvaule')
    for i in range(1,5):
        main(i)
    #main(1)
