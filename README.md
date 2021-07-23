# ExpertRNA

This algorithm usues a multi-branch reinforcement learning algorithm to predict secondary structures of peudoknot-free RNA sequences.

The current implementation uses RNAfold as the roll-out heuristic algorithm.

ENTRNA (with and without MFE) is used as expert.

### Environment
This requires the environment requirement listed in this repo:  
ENTRNA: https://github.com/sucongzhe/ENTRNA

We did some modifications on the original version of ENTRNA so the new version of ENTRNA is already inside the `util` folder. You don't need to download ENTRNA itself from the link above but you need the environments for ENTRNA from the link above. (There is a doc containing ENTRNA environment setting in the main folder)

We call ENTRNA expert with MFE MFE_ENTRNA and ENTRNA expert without MFE NFE_ENTRNA.

### Algorithms in this repository
There are three versions of this algorithom included in this repository: MFE only, NFE only and MFENFE_2X2.
MFE only is a 4 branch algorithm using only ENTRNA with MFE calculation as an expert.
NFE only is a 4 branch algorithm using only ENTRNA with NFE calculation as an expert.
MFENFE_2X2 is a 2 by 2 branches algorithm using ENTRNA with MFE and NFE as experts.

This respository contains two parts of work: (1) Running the alg iteself and produce an output file. (2) Exporting expert scores and free energy scores of the output file.

****************************
**For part 1:**

Running the whole algorithm requires three .py files:  
1. Batchrun_for_database_xxxxxxx.py  
2. Main_SeqAttRollOut_RNAfold_xxxxxxx.py
3. SeqAttRollOut_RNAfold_xxxxxxx.py

1 is the main python file to run. You may want to change the input file's directory and the name of output file within this script.

Input format:   
Example input format is as 5s_Acholeplasma-laidlawii-2_.ct file.

Output format:  
A csv file containing the following columns:

**Predictions**

1. 'Name', the name of the rna
2. 'Ori_Seq', the sequence of the RNA, notice that if your input RNA use T instead of U, we will automatically. transfer T into U
3. 'Actual_str', the actual structure of the input RNA, if any.
4. 'RNAfold_str', the RNAfold result structure of the input RNA.
5. 'Our_alg_str', the result structure by our alg of the input RNA.
6. 'RNAfold_distance', the distance percentage btw actual structure and structure by RNAfold
7. 'Our_distance', the distance percentage btw actual structure and structure by our alg

**The remaining columns are ENTRNA evaluations**

8. 'ent_3', sequence entropy with window size three, calculation details referred to Congzhe's paper
9. 'gc_percentage', the percentage of G and C in the input RNA seq.
10. 'ensemble_diversity', calculation details referred to Congzhe's paper
11. 'expected_accuracy', calculation details referred to Congzhe's paper
12. 'fe_per', calculation details referred to Congzhe's paper
13. 'Running_time(sec)', running time for this RNA in seconds


****************************
**For part 2:**

Exporting expert scores and free energy scores of the output file from part 1. This part is optional, you can make use of it if you want to know the scores of experts and free energies of the structures.

Input format:  
Just the csv file output from part 1.

Output format:  
The csv file from part 1 with the following columns added:  

**ENTRNA foldabilities**

14. 'Actual_Foldability_MFE', ENTRNA_MFE scores for the actual structure
15. 'Actual_Foldability_WOMFE', ENTRNA_NFE scores for the actual structure
16. 'Actual_abs_Foldability_difference', absolute difference of ENTRNA_MFE and ENTRNA_NFE scores for the actual structure

**ViennaRNA MFE scores**  

17. 'RNAfold_Foldability_MFE', ENTRNA_MFE scores for the structure by RNAfold
18. 'RNAfold_Foldability_WOMFE', ENTRNA_NFE scores for the structure by RNAfold
19. 'RNAfold_abs_Foldability_difference', absolute difference of ENTRNA_MFE and ENTRNA_NFE scores for the structure by RNAfold
20. 'Our_Foldability_MFE', ENTRNA_MFE scores for the structure by our alg
21. 'Our_Foldability_WOMFE', ENTRNA_NFE scores for the structure by our alg
22. 'Our_abs_Foldability_difference', absolute difference of ENTRNA_MFE and ENTRNA_NFE scores for the structure by our alg
23. 'Actual_FE', free energy (by. ViennaRNA) for the actual structure
24. 'RNAfold_FE', free energy (by. ViennaRNA) for the structure by RNAfold
25. 'Our_FE', free energy (by. ViennaRNA) for the structure by our alg




