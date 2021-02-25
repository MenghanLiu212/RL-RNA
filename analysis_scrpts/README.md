# rollout_plotter
These scripts provide analysis and visualization capabilities for ExpertRNA.  The ExpertRNA implementation in the parent folder produces as output a csv with the following columns:


| column number | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 | 21 | 22 | 23 | 24 | 25 |
|--------------|---|---|---|---|---|---|---|---|---|----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|
|**Data**|Number|Name|Sequence|Actual Structure|RNAfold Structure|ExpertRNA Structure|RNAfold Hamming Distance|ExpertRNA Hamming Distance|ent_3|GC Perentage|Ensemble Diversity|Expected Accuracy|fe_per|Running Time(sec)|Actual Foldability MFE|Actual Foldability NFE|Actual abs Foldability difference|RNAfold Foldability MFE|RNAfold Foldability NFE|RNAfold abs Foldability difference|ExpertRNA Foldability MFE|ExpertRNA Foldability NFE|ExpertRNA abs Foldability difference|Actual FE|RNAfold FE|ExpertRNA FE

## Script descriptions

**MCC.py**: Reads all .csv files in the current directory and calculates the [Matthews Correlation Coefficent](https://en.wikipedia.org/wiki/Matthews_correlation_coefficient) for every prediction compared with the actual structure.  This script also produces a large number of figures and statistics on the aggregate MCC data.  The user is encouraged to look through the plots and comment out those that they do not need.  The following files are produced:<br>
Three .fasta files for each line in the input .csv files: These files contain the structure name, the sequence, and each of the actual structure, the RNAfold prediction and the ExpertRNA prediction.  
Two .col files for each line in the input .csv files: These files contain a [Forna](http://rna.tbi.univie.ac.at/forna/) colormap marking which nucleotides were correctly predicted and which were incorrect.  
Five score files: One for each branch and one with the aggregate data (labeled 0) These have the name, and the MCC scores for the RNAfold prediction, and the ExpertRNA prediction and finally the difference between the two predictions.  This file is sorted from worst to best difference between ExpertRNA and RNAfold.

**forna_generator.py**: Reads all the .fasta and .col files produced by MCC.py in the present directory and spits out [Forna](http://rna.tbi.univie.ac.at/forna/) links for every structure.  These can be pasted in a web browser to view the actual, and predicted structures.  The predicted structures are annotated with red for incorrectly-predicted nuclotides and green for correct predictions.


## Usage
After running ExpertRNA, move the results file into its own directory.  Navigate to this directory and run
```shell
python /path/to/RL-RNA/analysis_scripts/MCC.py
```

This will produce a number of PNG-format overview graphs, FASTA structure files and Forna colormap files.

Copy the files applicable to strucrures that you are particularly interested in (use the scores_X.dat files to select particularly good or bad structures) to a fresh directory and run
```shell
python /path/to/RL-RNA/analysis_scripts/forna_generator.py
```

This will produce a forna link for each structure in the folder. Copy-paste these into your web browser to view the RNA structure.