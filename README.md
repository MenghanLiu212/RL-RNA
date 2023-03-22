# ExpertRNA

This algorithm usues a multi-branch reinforcement learning algorithm to predict secondary structures of peudoknot-free RNA sequences.

The current implementation uses RNAfold as the roll-out heuristic algorithm.

ENTRNA (with and without MFE) is used as the expert.

### Environment
The Anaconda environment which was used to test and run ExpertRNA can be found in the expertRNA_environment.txt file in this repository.  Please note that the ViennaRNA package (including RNAfold) must be compiled from source on Python versions >= 3.8.  We use 3.7 for testing as ViennaRNA can be directly installed from Bioconda for that version.

<ins> Environment summary </ins>
* Python 3.7  
* numpy
* scikit-learn
* pandas
* viennarna

We made some modifications to the original version of ENTRNA, so the new version of ENTRNA is already inside the `util` folder. You don't need to download ENTRNA itself from the link above but you need the environments for ENTRNA from the link above. (There is a doc containing ENTRNA environment setting in the main folder)

We call ENTRNA expert with MFE included in the parameterization MFE_ENTRNA, and ENTRNA expert without MFE in the parameterization NFE_ENTRNA.

Training dataset:
 - ENTRNA: 1595 sequences; from the RNASTRAND database: RL-RNA/util/RNASTRAND_pseudoknot_free_feature.csv
    
Testing dataset:
> Each dataset was culled to remove structures with greater than 80% sequence identity withstructures in the ENTRNA training dataset using CD-HIT-EST-2D (Huang et al. 2010).
 - ExpertRNA: 
    - Mathews data: 1559 sequences; RL-RNA/testing_data/mathews (the folder contains xxx.dbn files each of which contains a sequence-structure pair)
    - rfam data: 147 sequences; RL-RNA/testing_data/rfam (the folder contains xxx.dbn files each of which contains a sequence-stricture pair)

## References
Huang, Y., Niu, B., Gao, Y., Fu, L., & Li, W. (2010). CD-HIT Suite: a web server for clustering and comparing biological sequences. Bioinformatics, 26(5), 680-682.

Liu, M., Poppleton, E., Pedrielli, G., Šulc, P., & Bertsekas, D. P. (2022). Expertrna: A new framework for rna secondary structure prediction. INFORMS Journal on Computing, 34(5), 2464-2484.

****************************

### Implementation
Running the whole algorithm requires three .py files:  
1. ExpertRNA.py  
2. ExpertRNA_main.py
3. ExpertRNA_toolbox.py

1 is the main python file to run. You may want to change the input file's directory and the name of output file within this script. 

There are two modes of this software: test mode and prediction mode. For both mode, the input format can be either .dbn (for testing) or .fasta (for prediction). The output of prediction mode will be a .dbn file with predcited strucutres in dot-bracket format and corresponding expert scores for each branch. The output of test mode is a csv containing the following columns: 

**Predictions**

1. 'Name', the name of the rna
2. 'Ori_Seq', the sequence of the RNA, notice that if your input RNA use T instead of U, we will automatically convert T into U
3. 'Actual_str', the actual structure of the input RNA, if known.
4. 'RNAfold_str', the RNAfold result structure of the input RNA.
5. 'Our_alg_str', the result structure by ExpertRNA of the input RNA.
6. 'RNAfold_distance', the distance percentage between the actual structure and structure by RNAfold
7. 'Our_distance', the distance percentage between actual structure and structure by ExpertRNA

**The remaining columns are ENTRNA evaluations**

8. 'ent_3', sequence entropy with window size three, calculation details referred to Congzhe's paper
9. 'gc_percentage', the percentage of G and C in the input RNA seq.
10. 'ensemble_diversity', measures the expected distance between the target secondary structure and all the other secondary structure, calculation details referred to Congzhe's paper.
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
19. 'Actual_FE', free energy (from ViennaRNA) for the actual structure
20. 'RNAfold_FE', free energy (from ViennaRNA) for the structure by RNAfold
21. 'Our_FE', free energy (from ViennaRNA) for the structure by ExpertRNA


### Commandline instructions

This software is run via commandline. The options are as follow:

1. `input_data`: The path to a directory containing input sequences. Each sequence should be a separate fasta or dbn file.
2. `-e`: The Expert used to evaluate partial solutions, and branches to maintain for each expert. Expert should choose from 'ENTRNA_MFE', 'ENTRNA_NFE'. Should be in the form of `expert_name branch_num`.  You may call `-e` multiple times to run multiple branches with different experts.
3. `-f`: The folder to complete partial solutions. Currently the only option is 'RNAfold', but more will be added. Type should be 'nonspecific'. Should in the form of `folder_name type`.
4. `-o`: Name of file to write results out to in .csv or .dbn format.
5. `-m`: The minimum basepair distance, default is 4.
6. `-t`: Testing flag.  If set, the algorithm expects inputs to be in .dbn format (each file has 3 lines, a fasta-style header, the sequence, and a corresponding correct structure) and the output will be in the .csv format explained in [Implementation](#Implementation)

Example:
For input data in DatasetPath with expert as 'ENTRNA_MFE' and folder as 'RNAfold', with minimum basepair distance as 4, in test mode, the command line invocation is:
```bash
python expertRNA.py input_data -e ENTRNA_MFE 4 -f RNAfold nonspecific -m 4 -t
```
