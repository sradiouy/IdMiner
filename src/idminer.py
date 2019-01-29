import os
import logging
import time
import FetchArticles
import TermsByAbstracts
import CleanTerms
import IdMinerReports

#In order to run....

#Need to load form dash app:

    # Input_file: txt or fasta file.
    # Format_file: selection of text or fasta (radio button in app)
    # Run_name: name of output files (text box in app) 

#Example of run:

Input_file = os.path.join(os.path.dirname(__file__), "../data/Test_Fetch_articles/Guille_genes.txt")
Format_file = "text"
Run_name = "Guille_run"
start_time = time.time()
logging.basicConfig(format='%(asctime)s - IdMiner - %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p',filename= Run_name + '.log',level=logging.DEBUG)
FetchArticles.ids_by_gene(Run_name,Input_file,Format_file,30,30)
elapsed_time = time.time() - start_time
logging.info("Fetch article -> Duration %i seconds" %(elapsed_time))
start_time = time.time()
gene_pubmed = TermsByAbstracts.gene_articles_dict(Run_name + ".tsv")
abstractdict = TermsByAbstracts.generate_abstracts_dict(gene_pubmed)
elapsed_time = time.time() - start_time
logging.info("Abstract by Term -> Duration %i seconds" %(elapsed_time))
start_time = time.time()
worddict = CleanTerms.cleanwords(abstractdict)
dictgram = CleanTerms.dict_gram(worddict,mostcommon=5000,minfreq=2)
elapsed_time = time.time() - start_time
logging.info("Clean Terms -> Duration %i seconds" %(elapsed_time))
start_time = time.time()
termdict,geneterms, geneNterms = IdMinerReports.get_gene_term_dicts(gene_pubmed,worddict,dictgram)
#IdMinerReports.create_co_occurrence_matrix(gene_pubmed,dictgram,geneNterms,Run_name)
IdMinerReports.create_terms_info_dataframe(termdict,Run_name)
IdMinerReports.create_gene_artciles_dataframe(geneterms,Run_name)
elapsed_time = time.time() - start_time
logging.info("Reports -> Duration %i seconds" %(elapsed_time))
