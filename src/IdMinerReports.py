from collections import defaultdict
import pandas as pd
from src.IdfZipfScore import *
from src.TermsByAbstracts import check_df_ids_by_gene, gene_articles_dict


def get_gene_term_dicts(genepmids,worddict,dictgram):
    """Returns a collection of dicts. With infomration of the term and its association with genes. 
    
    Arguments:
        genepmids {[dict]} -- Dictionary of {genes id: pubmed ids}
        worddict {[type]} -- Dictionary of articles, where the key is the pubmed id and the value is a list of cleaned tokens.
        dictgram {[type]} -- Dictionary of terms (key). Values are [freq,# articles, pubmed ids]
    
    Returns:
        [dict] -- termdict: dictionary of terms. Contain information of ZIPF, IDF, Freq, #Articles, Publications. geneterms: dictionary of terms and genes. Contain publications were the term appear associated with the gene. geneNterms: the count of the previous dictionary.
    """

    termdict = defaultdict(dict) # term information dict 
    geneterms = defaultdict(dict) # term by gene, pubmedid
    geneNterms = defaultdict(dict) # term by gene, number of pubmedid
    zipf = zipf_score(dictgram) #dictionary of zipf scores
    idf = idf_dict(worddict,dictgram) #dictionary of idf scores
    gramlist = dictgram.keys() #get list of terms
    for term in gramlist:
        for gene in genepmids:
            geneterms[term][gene] = ",".join([x for x in dictgram[term][2].split(",") if x in genepmids[gene].split(",")])
            if geneterms[term][gene] != "":
                geneNterms[term][gene] = geneterms[term][gene].count(",") +1 
            else:
                geneNterms[term][gene] = 0
        termdict[term]["Articles"] = dictgram[term][1] #Number of articles where the term appear.
        termdict[term]["Genes"] = sum([1 for gene in geneterms[term] if geneterms[term][gene] != ""]) #Number of genes associated with the term 
        termdict[term]["Freq_Term"] = int(dictgram[term][0]) #Freq of term in corpus
        termdict[term]["ZIPF_Score"] =  zipf[term] #zipf score (general freq of the word)
        idf_score =round(float(idf[term]),3)
        termdict[term]["IDF_Score"] = idf_score #idf score
        termdict[term]["Publications"] = dictgram[term][2] # publications
    return  termdict,geneterms, geneNterms   

def create_co_occurrence_matrix(genepmids,dictgram,geneNterms,run_name):
    """ Create a co-occurrence matrix of genes and terms. Each term has a number of articles associated with each gene.
    
    Arguments:
        genepmids {[dict]} -- Dictionary of {genes id: pubmed ids}
        dictgram {[dict]} -- Dictionary of terms (key). Values are [freq,# articles, pubmed ids]
        geneNterms {[dict]} -- Dictionary of terms and genes. Contain number of publications were the term appear associated with the gene.
        run_name {[string]} -- Run name for output files
    """

    dfco = pd.DataFrame(columns = list(genepmids.keys()), index = dictgram.keys())
    dfco[:] = int(0)
    for gene in list(genepmids.keys()):
        for term in dictgram:
            dfco[gene][term] = geneNterms[term][gene]
    dfco.to_csv("matrix-"+run_name+".csv",index=True,header=True,sep=",")

def create_terms_info_dataframe(termdict,run_name):
    """ Create a dataframe of infromation about terms.  
    
    Arguments:
        termdict {[dict]} -- termdict: dictionary of terms. Contain information of ZIPF, IDF, Freq, #Articles, Publications.
        run_name {[string]} -- Run name for output files
    """

    dfterms = pd.DataFrame(termdict).T #transpose
    dfterms["Terms"] = dfterms.index.tolist() #index to column name
    dfterms = dfterms[["Terms","Genes","ZIPF_Score","Articles","Freq_Term","Publications"]] #change order
    dfterms.to_csv(run_name+"_IDMiner-Terms.csv",index=False,header=True,sep=",") # to csv

def create_gene_artciles_dataframe(geneterms,run_name):
    """ Create a datframe of genes and terms. Each term has articles associated with each gene.
    
    Arguments:
        geneterms {[dict]} -- Dictionary of terms and genes. Contain publications were the term appear associated with the gene.
        run_name {[string]} -- Run name for output files
    """

    dfgenesbyarticles =pd.DataFrame(geneterms).T
    dfgenesbyarticles["Terms"] = dfgenesbyarticles.index.tolist()
    dfgenesbyarticles.to_csv(run_name+"_IDMiner-Genes.csv",index=False,header=True,sep=",")

#Crear un solo csv de los tres. 

# df.to_sql