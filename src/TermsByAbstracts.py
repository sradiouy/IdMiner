import pandas as pd
import backoff
import logging
import requests
from Bio import Entrez 
from src.FetchArticles import generate_random_emails
from src.HelperGrams import generate_stop_grams, get_keepterms


def check_df_ids_by_gene(articles):
    """ Control input file to be processed by IdMiner
    
    Arguments:
        input_file {[file]} -- tsv file containing gene ids and pubmed ids associated with it
    
    Raises:
        FileNotFoundError -- Raise an error when file does not exist
        ValueError -- Raise an error when the number of column's is different to 2
        ValueError -- Raise an error when the names of columns are not [id,article]
    
    Returns:
        [dataframe] -- Panda dataframe (tsv format). Two columns [genes, pubmed ids]
    """

    # try:
    #     df = pd.read_csv(input_file,sep="\t",header=0)
    #     df.dropna(inplace=True) # Remove empty genes
    # except FileNotFoundError:
    #     mssg = "There was not file named %s" %(input_file)
    #     logging.error(mssg)
    #     raise FileNotFoundError(mssg)
    df = pd.DataFrame(articles,columns=["id","article"])
    # if len(df.columns) != 2:
    #     mssg = "%s has wrong number of columns. Must have to be named (%s - %s), and must be tsv format. Your %s has %i column's" % (input_file,"id","article",input_file,len(df.columns))
    #     logging.error(mssg)
    #     raise ValueError(mssg)
    # elif df.columns[0] != "id" or  df.columns[1] != "article":
    #     mssg = "%s has wrong name of columns. Must have to named (%s - %s), and must be tsv format." % (input_file,"id","article")
    #     logging.error(mssg)
    #     raise ValueError(mssg)
    return df 

def gene_articles_dict(articles):
    """ Take a tsv file of genes and pubmed ids and transform it into a dictionary.
    
    Arguments:
        input_file {[file]} -- tsv file containing gene ids and pubmed ids associated with it
        part {[string]} -- search in abstract or full-text body
    
    Returns:
        [dict] -- Dictionary of {genes id: pubmed ids}
    """

    df = check_df_ids_by_gene(articles)
    genepmids = {}
    for r in df.iterrows():
        gene = r[1].id
        pmids = str(r[1].article)
        genepmids[gene] = pmids
    return genepmids

@backoff.on_exception(backoff.expo,requests.exceptions.RequestException,max_time=120)
def fetch_abstarcts(gene,pubmed,abstractdict):
    """Recover abstracts from pubmed ids through the use of the NCBI efetch function of the Entrez Biopython module.
    
    Arguments:
        gene {[string]} -- [description]
        pubmed {[string]} -- [description]
        abstractdict {[type]} -- [description]
    
    Returns:
        [dict] -- Dictionary of articles, where the key is the pubmed id and the value is the abstract
    """
    Entrez.email = generate_random_emails()
    logging.info("Fetching abstracts for: %s" % (gene))
    nabstracts = 0
    try:
        handle = Entrez.efetch(db="pubmed", id=pubmed,rettype="xml", retmode="text",retmax=100000) # From NCBI obtain information of pubmed articles
        records = Entrez.read(handle) # retain records 
        passgene = False
    except:
        passgene = True
        logging.warning("There was an error with %s can't fetch NCBI articles. Skipping %i..." %(gene,pubmed.count(",")+1))
    if not passgene:
        for pubmed_article in records['PubmedArticle']: 
            try:
                ids = pubmed_article['MedlineCitation']['PMID'].rstrip()
                abstract = ".".join(map(str,pubmed_article['MedlineCitation']['Article']['Abstract']['AbstractText']))
                abstractdict[ids] =abstract         
            except:
                logging.warning("There isn't abstract for pubmed article %s. Skipping..." %(ids))
                continue
            nabstracts += 1
        logging.info(" -- we can retrieve %i articles" %(nabstracts))
    return abstractdict

def generate_abstracts_dict(gene_pubmed):
    """ Generates a dictionary of articles, where the key is the pubmed id and the value is the abstract
    
    Arguments:
        gene_pubmed {[dict]} -- Dictionary of {genes id: pubmed ids}
    
    
    Returns:
        [dict] -- Dictionary of articles, where the key is the pubmed id and the value is the abstract 
    """

    abstractdict = {}
    for gene,pubmed in gene_pubmed.items():
        abstractdict = fetch_abstarcts(gene,pubmed,abstractdict)
    logging.info("There was in total %i unique articles for %i genes" %(len(abstractdict),len(gene_pubmed)))
    return  abstractdict


