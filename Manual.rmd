# Usage: 

1. Setup 

IdMiner expects a list of gene IDs (as a text file with .txt extension) or a multifasta file (.fasta extension).  

The user may establish specific terms to be kept (even if they are not among the 5000 most frequent terms) or to be excluded (even if they are common English words). The frequency (i.e., high frequency means that the word analyzed includes words that are frequent in English, like human or cancer) and the number of terms kept can be modified. Finally, the Coverage and Identity parameters are used during the blast step to define homolog genes 

2. Explore the results 

Upon completion of the run, results will be stored in two files located in IdMiner/project/Results. You may move these files to other locations. Files will be named “XXX_IDMiner-Genes.csv” and “XXX_IDMiner-Terms.csv” where XXX corresponds to the name of the GI list file that was analyzed 

To browse the results, go to either the "Term Exploration" or the "Gene Exploration" tabs.  In the "Term Exploration" tab, open both the XXX_IDMiner-Genes.csv and the XXX_IDMiner-Terms.csv file. Finally, in the "Gene Exploration" tab, open only the XXX_IDMiner-Genes.csv file. 

 
 
