# ExpertRNA

This algorithm usues a multi-branch reinforcement learning algorithm to predict secondary structures of peudoknot-free RNA sequences.

The current implementation uses RNAfold as the roll-out heuristic algorithm.

ENTRNA (with and without MFE) is used as expert.

### Environment
This requires the environment requirement listed in this repo:  
ENTRNA: https://github.com/sucongzhe/ENTRNA

We did some modifications on the original version of ENTRNA so the new version of ENTRNA is already inside the `util` folder. You don't need to download ENTRNA itself from the link above but you need the environments for ENTRNA from the link above. (There is a doc containing ENTRNA environment setting in the main folder)

We call ENTRNA expert with MFE MFE_ENTRNA and ENTRNA expert without MFE NFE_ENTRNA.


****************************

### Implementation
Running the whole algorithm requires three .py files:  
1. ExpertRNA.py  
2. ExpertRNA_main.py
3. ExpertRNA_toolbox.py

1 is the main python file to run. You may want to change the input file's directory and the name of output file within this script. 

There are two modes of this software: test mode and prediction mode. For both mode, the input format can be both .fasta and .dbn files. The output of prediction mode will be a .dbn file with predcited strucutres in dot-bracket format and corresponding expert scores for each branch. And the output of test mode is a csv containing the following columns: 

**Predictions**

1. 'Name', the name of the rna
2. 'Ori_Seq', the sequence of the RNA, notice that if your input RNA use T instead of U, we will automatically. transfer T into U
3. 'Actual_str', the actual structure of the input RNA, if any.
4. 'RNAfold_str', the RNAfold result structure of the input RNA.
5. 'Our_alg_str', the result structure by ExpertRNA of the input RNA.
6. 'RNAfold_distance', the distance percentage btw actual structure and structure by RNAfold
7. 'Our_distance', the distance percentage btw actual structure and structure by ExpertRNA

**The remaining columns are ENTRNA evaluations**

8. 'ent_3', sequence entropy with window size three, calculation details referred to Congzhe's paper
9. 'gc_percentage', the percentage of G and C in the input RNA seq.
10. 'ensemble_diversity', measures the expected distance between the target secondary structure and all the other secondary structure, calculation details referred to Congzhe's paper
11. 'expected_accuracy', measures the expected number of bases that are in correct base pairing status, calculation details referred to Congzhe's paper
12. 'Running_time(sec)', running time for this RNA in seconds 

**ENTRNA foldabilities**

13. 'Actual_Foldability_MFE', ENTRNA_MFE scores for the actual structure
14. 'Actual_Foldability_NFE', ENTRNA_NFE scores for the actual structure

**ViennaRNA MFE scores**  

15. 'RNAfold_Foldability_MFE', ENTRNA_MFE scores for the structure by RNAfold
16. 'RNAfold_Foldability_NFE', ENTRNA_NFE scores for the structure by RNAfold
17. 'Our_Foldability_MFE', ENTRNA_MFE scores for the structure by ExpertRNA
18. 'Our_Foldability_NFE', ENTRNA_NFE scores for the structure by ExpertRNA
19. 'Actual_FE', free energy (by. ViennaRNA) for the actual structure
20. 'RNAfold_FE', free energy (by. ViennaRNA) for the structure by RNAfold
21. 'Our_FE', free energy (by. ViennaRNA) for the structure by ExpertRNA


### Implementation commandline instructions

This software is run by commandline. The options are as follow:

1. 'input_data': The path to a directory containing input sequences. Each sequence should be a separate fasta file.
2. '-e': The Expert used to evaluate partial solutions, and branches to maintain for each expert. Expert should choose from 'ENTRNA_MFE', 'ENTRNA_NFE'. Should in the form of [expert name branch num]. And input should be call the -fd multiple times for each input.
3. '-f': The folder to complete partial solutions. Should choose from 'RNAfold' in this version. Type should be nonspecific. Should in the form of [folder name type]. And input should be call the -fd multiple times for each input.
4. '-o': Name of file to write results out to in .csv format.
5. '-m': The minimum basepair distance, default is 4.
6. '-t': If set the algorithm expects a true structure for each input on the line following the sequence and produces a comparison between the prediction and the real score.

Example:
For input data in DatasetPath with expert as 'ENTRNA_MFE' and folder as 'RNAfold', with minimum basepair distance as 4, in test mode:
Type in: $ python ExpertRNA_v3.py DatasetPath -e 'ENTRNA_MFE' 4 -f 'RNAfold' 'nonspecific' -m 4 -t
