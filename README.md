# RL-RNA

This algorithm is using a multi-branch reinforcement learning algorithm to predict secondary structure of peudoknot-free RNA sequence.
RNAfold is used as the roll-out heuristic algorithm.
ENTRNA (with and without MFE) is used as expert.
We call ENTRNA expert with MFE MFE_ENTRNA and ENTRNA expert without MFE NFE_ENTRNA.

There are three versions of this algorithom included in this repository: MFE only, NFE only and MFENFE_2X2.
MFE only is a 4 branch algorithm using only ENTRNA with MFE calculation as an expert.
NFE only is a 4 branch algorithm using only ENTRNA with NFE calculation as an expert.
MFENFE_2X2 is a 2 by 2 branches algorithm using ENTRNA with MFE and without MFE (we call it NFE) as experts.

Running the whole algorithm requires three .py files:
(1) Batchrun_for_database_xxxxxxx.py
(2) Main_SeqAttRollOut_RNAfold_xxxxxxx.py
(3) SeqAttRollOut_RNAfold_xxxxxxx.py
And the (1) is the main python file to run. You may want to change the input file's directory and the name of output file.


Imput format: 
Example input format is as 5s_Acholeplasma-laidlawii-2_.ct file.

