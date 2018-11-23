#!/usr/bin/env python3.6
import os
import sys
from pathlib import Path
import subprocess
from bs4 import BeautifulSoup as BS
import urllib.request
import itertools
from wordfreq import zipf_frequency
from lxml import etree as ET
import pandas as pd
import urllib
import numpy as np
import argparse
import re
from Bio import Entrez 
import time
import csv
import backoff
import importlib
import nltk
import nltk.corpus
import numpy as np
import math
from textblob import TextBlob as tb
from PIL import Image
from os import path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import random
from wordcloud import WordCloud, STOPWORDS, ImageColorGenerator
import becas
import json
from gooey import Gooey, GooeyParser
from collections import Counter
from collections import defaultdict
from nltk.tokenize import RegexpTokenizer
from nltk.stem.snowball import SnowballStemmer
import string



# running = True
# nltk.download("stopwords", quiet=True)
# nltk.download("punkt", quiet=True)
# nltk.download("wordnet", quiet=True)

@Gooey(advanced=True,
    optional_cols=2,
    program_name="Gene Cloud",
    tabbed_groups=True,
    dump_build_config=True,
    show_success_modal=False,
    default_size=(1024, 768))

def window():
    parser = GooeyParser(description="An enhancer for PaperBlast")
    subs = parser.add_subparsers(help='commands', dest='command')
    parser_one = subs.add_parser('Module 1', prog="PaperBlast")
    parser_one.add_argument('genes',default="id1\nid2", widget='Textarea',help="Uniprot IDs of your genes of interest. Separate by new lines",
                        metavar='Uniprot IDs',gooey_options={
                            'height': 10000})
    parser_one.add_argument("output",action="store",default="output.txt",widget='TextField',
                        metavar='Output File',help="Name of the output file")
    parser_two = subs.add_parser('Module 2', prog="Enhancer")
    parser_two.add_argument('fasta',metavar='multifasta file',action='store',widget='FileChooser',help="multifasta file in fasta format")
    parser_two.add_argument('pubmed',metavar='Pubmed IDs files',action='store',widget='FileChooser',help="Output file of Woody 1")
    parser_two.add_argument('basename', default="basename", widget="TextField",help="basename of the output files.",metavar='Basename')
    parser_two.add_argument('--terms', default="regex1\nregex2", widget="Textarea",help="Regex terms that the abstracts must have",
                        metavar='Terms',gooey_options={
                            'height': 10000})
    listimg = os.listdir(os.path.join(os.path.dirname(__file__), './pics'))
    parser_two.add_argument('--img',
                        metavar="Wordcloud Image",
                        help = "Choose one",
                        nargs="+",
                        default=['square'],
                        choices=sorted(listimg),
                        widget='Listbox',
                        gooey_options={
                            'height': 1000,
                            'heading_color': 'blue',
                            'text_color': 'red',
                        }
                        )
    parser_two.add_argument('--stopwords', default="word1\nword2",widget='Textarea',help="Stopword to clean your results",
                        metavar='Stopwords')
    parser_two.add_argument('--pubtator', default="", widget="TextField",help="output file name,  if empty not run Pubtator",metavar='Pubtator results')
    parser_two.add_argument('--project', default="", widget='TextField',help="Tagtog project name",metavar="Project name")
    parser_two.add_argument('--username', default="", widget='TextField',help="Tagtog username",
                        metavar='Username')
    parser_two.add_argument('--password', default="", widget='PasswordField',help="Tagtog password",
                        metavar='Password')
    args = parser.parse_args()
    return args

#PaperBlast

# El primer paso del programa es buscar los genes en el paperBlast. Para eso utilizamos dos formas distintas. 
# En esta primera funcion se busca a partir del identificador. Se tiene en cuenta varios intentos de busqueda.


#### Modulo 1 #####

@backoff.on_exception(backoff.expo,urllib.error.URLError,max_value=32)
def get_text(gene):
    """
    Esta funcon permite ejecutar una busqueda de paperblast y obtener el resultado en formato html.
    input: gene: nombre de un gen.
    output: list: lista en formato en html, donde cada elemento esta conformado por un hit.
    """
    try:
        print("Searching %s in paperblast!" % (gene))
        primarylist = []
        url = 'http://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi?query=' + gene
        response = urllib.request.urlopen(url)
        data = response.read()      # a `bytes` object
        text = data.decode('utf-8')
        list = []
        for x in str(text).split("\n"):
            if "showAlign" in x:
                list.append(x)
        return list
    except:
        return False

# El primer paso del programa es buscar los genes en el paperBlast. Para eso utilizamos dos formas distintas.
# En esta segunda funcion se busca a partir del identificador. Se tiene en cuenta varios intentos de busqueda. 

@backoff.on_exception(backoff.expo,urllib.error.URLError,max_value=32)
def get_text_from_fasta(seq):
    try:
        print("Searching gene in paperblast!")
        primarylist = []
        url = 'http://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi?query=' + seq
        response = urllib.request.urlopen(url)
        data = response.read()      # a `bytes` object
        text = data.decode('utf-8')
        list = []
        for x in str(text).split("\n"):
            if "showAlign" in x:
                list.append(x)
        return list
    except:
        return False

#Para cada gen se busca obtener la informacion del blast. Esta incluye infromacion de como fue el hit vs los homologos (id,coverage,bitscore y evalue) y en que articulo fue nombrado cada gen homologo.

def get_info_hit(results_list,pcov,pident):
    """
    La siguiente funcion me permite a partir del resultado de una ejecucion de paperBlast obtener un resumen de esa busqueda. Indicando para cada hit de nuestro gen de interes, los parametros de blast y los articulos deseados.
    -input: result_list: html resultado de una ejecucion de blast donde cada hit esta en un elemento de la lista. pcov: porcentaje de cobertura minimo para que el hit sea aceptado. pident: porcentaje de identidad minimo para que el hit sea aceptado.
    -Out: hit es un diccionario donde tenemos toda la informacion parseada del html. Para cada gen.
    """
    regex = "title=.*?(.+?)</a>"
    hit = {}
    count = 0
    for x in results_list:
        results = re.findall(regex,x)
        for n,y in enumerate(results):
            if n == 0:
                db = y.split(">")[0].split(" ")[0].replace('"',"")
                gene = y.split(">")[1].split(" ")[-1]
            else:
                if "% identity, " in y:
                    values = y
                    evalue,score = values[values.find("(")+1:values.find(")")].split(",")
                    evalue = evalue.replace("E = ","")
                    score = score.replace(" bits","").strip()
                    identity,coverage = values.split(">")[-1].split(",")
                    identity = identity.split("%")[0].strip()
                    coverage = coverage.split("%")[0].strip()
                    continue
        if pcov <= int(coverage) and pident <= int(identity):
            soup = BS(x,"lxml")
            primarylist = []
            for a in soup.findAll('a',href=True):
                if re.findall('pmc|pubmed', a['href']):
                    primarylist.append(a['href'])
            if primarylist != []:
                count += 1
                key = "Hit_" + str(count)
                hit[key] = db,gene,evalue,score,identity,coverage,primarylist
        else:
            pass
    return hit

@backoff.on_exception(backoff.expo,urllib.error.URLError,max_value=2)
def convert_pmc(pmc):
    """
    Permite convertir los ids de pmc a pubmed utuilizando una api de ncbi. Como input una lista de ids de pmc.
    input: pmc: lista de id de pmc.
    output: lista de ids convertidos a pubmed id.
    """
    #print("Converting PMC ids to PMID!")
    pmclist = []
    for x in pmc:
        try:
            url = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?tool=genomepaper&email=anonymus@gmail.com&ids=" + x
            response = urllib.request.urlopen(url)
            data = response.read() 
            soup = BS(data,"lxml")
            for a in soup.findAll('record',pmid=True):
                pmclist.append(a['pmid'])
        except:
            pass
    return pmclist

def get_list(hits):
    primarylist = []
    for hit in hits:
        primarylist += hits[hit][-1]
    listpmc = []
    listpubmed = []
    for x in primarylist:
        pmc = [s.rstrip() for s in x.split("/") if "PMC" in s]
        if len(pmc) != 0:
            listpmc.append(pmc[0])
        else:
            if "pubmed" in x:
                listpubmed += (x.split("pubmed/")[-1].split(","))
    listpmc = list(filter(None,listpmc))
    listpubmed = list(filter(None,listpubmed))
    listpmc = list(set(listpmc))
    listpubmed = list(set(listpubmed))
    listpmc = ",".join(listpmc).split(",")
    pmc = convert_pmc(listpmc)
    listpubmed = ",".join(listpubmed).split(",")
    pmids = pmc + listpubmed
    pmids = list(set(pmids))
    pmids = list(filter(None,pmids))
    return ",".join(pmids)

def get_homology(hit,gene):
    homology = []
    for key in hit.keys():
        text_list = [gene,hit[key][1],hit[key][0],hit[key][2],hit[key][3],hit[key][4],hit[key][5],",".join((map(str,hit[key][6])))]
        homology.append(text_list)
    return homology

def create_paper_blast_table(genes):
    articles = []
    homology = []
    for x in genes.split("\n"):
        all_hits = get_text(x)
        if all_hits == []:
            pass
        else:
            hit_info = get_info_hit(all_hits,30,30)
            tmphomology = get_homology(hit_info,x)
            homology += tmphomology
            pmids = get_list(hit_info)
            tmpart = [x] + [pmids]
            articles += [tmpart]
    dfhomology = pd.DataFrame(homology,columns = ["gene","homologue","dataBase","evalue","bitscore","identity","coverage","ids"])
    dfarticles = pd.DataFrame(articles,columns=["id","article"])
    dfhomology.to_csv("Homology.tab",sep="\t",header=True,index=False)
    dfarticles.to_csv("Articles.tab",sep="\t",header=True,index=False)
    return None

################################
####          Modulo 2     #####    
################################

def ids_by_gene(df):
    """
    Function used to show the articles that each gene has in the paperblast search.
    input: df: dataframe of articles (two cols, id -gene- and articles -pubmedid-)
    output: genepmids: dictionary of gene (key) and articles (value).
    """
    genepmids = {}
    for index,row in df.iterrows():
        gene = row["id"]
        pmids = row["article"]
        genepmids[gene] = pmids
    return genepmids

@backoff.on_exception(backoff.expo,urllib.error.URLError,max_value=2)
def get_abstarcts(pmids,gene,abstractdict,email="anmus@gmail.com"):
    """
    Funcion que permite tomar una lista de articulos (pubmedids) y devuelve una lista de abstarcts.
    input: pmids: lista de articulos. gene: nombre del gen analizado.
    output: abstarcts: lista de abstracts.
    """
    Entrez.email = email
    print("Analyzing %s" % (gene))
    handle = Entrez.efetch(db="pubmed", id=','.join(map(str, pmids.split(","))),rettype="xml", retmode="text")
    records = Entrez.read(handle)
    for pubmed_article in records['PubmedArticle']:
        try:
            ids = pubmed_article['MedlineCitation']['PMID'].rstrip()
            abstract = ".".join(map(str,pubmed_article['MedlineCitation']['Article']['Abstract']['AbstractText']))
            if abstract == " " and ids != " ":
                print(ids,abstract)
            abstractdict[ids] =abstract          
        except:
            print(ids, "abstract not found!")
            continue
    return abstractdict


def genes_list_abstarcts(genepmids):
    """
    Function to obtain the one and bigrams associated with a specific gene. 
    input: pmids: list of pubmed id. gene: gene who is associated with the pmids.
    output: genedict: dictionary of genes and his associated one and bigrams.
    """
    abstractdict = {}
    for gene,pmids in genepmids.items():
        abstractdict = get_abstarcts(pmids,gene,abstractdict,"anmus@gmail.com")
    return  abstractdict


def get_words(abstractdict):
    worddict = {}
    removepunct = string.punctuation.replace("-","").replace("_","")
    table = str.maketrans('', '', removepunct)
    wordsbyabstracts = []
    for article in abstractdict: 
        worddict[article] = ([*map(str.lower, nltk.word_tokenize(abstractdict[article].translate(table)))])
    return worddict


#biology_common_words 
#"https://raw.githubusercontent.com/glutanimate/wordlist-medicalterms-en/master/wordlist.txt"

#Calcule la frecuencia, y esas son las que son mayores a 3.4 y significativas. 

def get_keepterms():
    keepterms = ['head','growth', 'oil','age','pain','air','sleep','dead', 'mother','human', 'earth', 'body', 'hand','network', 'island','sea','face','baby','heart', 'population', 'river', 'wall', 'pressure', 'weight','food','blood','cat','dog','nail', 'defensive','potato','pepper','bird', 'grass','marijuana','translation','oxygen', 'california', 'indian', 'atlanta', 'pig', 'intelligence', 'silence','foot','nonsense','shell', 'rabbit', 'ghost', 'elephant','ethnic', 'stress', 'clock','compound', 'egyptian','genius', 'manchester', 'breathe', 'oral', 'hostile','vacuum', 'fox', 'tension','humor', 'beer', 'metal','circuit','nest','philadelphia','mineral', 'thermal','snow','widow', 'bug', 'cognitive','blind','blast', 'patient','japanese','inflation','nutrition','signature', 'breed','buffalo', 'cope','trunk', 'giant', 'rat', 'belt','bread', 'pole','protection', 'tunnel', 'hmm', 'attract', 'brain','skull', 'solid', 'depression','taste','brussels','racial', 'survival','deer', 'queensland', 'cave','toronto', 'transfer', 'smell','skin', 'shock', 'brazilian','dutch', 'height', 'essence', 'gdp','sweet','salt','arizona','geneva','wire','columbia', 'cancer','environment', 'lip', 'asian', 'burst', 'throat','mountain', 'micro', 'invasion','diabetes', 'boundary','zone','kenya','rhythm','shield','bone','pregnancy', 'beauty', 'todd','ring','honey','toe', 'memory','arc', 'flower','nerve', 'vision','eagle', 'wind','italian','transition','dominant', 'bite', 'flash','orientation','horse', 'virus','repair','spider', 'tropical', 'berlin', 'crown','sexual', 'kidney', 'warm', 'resistance','canal','integrity', 'chile','progressive','pregnant','autumn','hook', 'hole', 'disk', 'sheep','breast', 'root','childhood','coat', 'viral','arabia', 'behavior', 'victoria','vienna','deep', 'communicate','arctic', 'shark', 'fitness', 'urban', 'african','grain','holland','glasgow', 'rice','wheat', 'kennedy', 'matrix','female', 'britain', 'flame', 'timber', 'wing','pine', 'motivation','moon', 'cure', 'disease', 'olive','stomach', 'bee','breath','fat','electric', 'bucket','gulf', 'consciousness', 'argentina','cow','panic','cycle','visual','habitat','gap','dublin','diet','delhi','teeth','bell','male','plant','leaf','facial','newcastle','lamb','marine','sydney','conscience', 'mirror','apple', 'killer', 'duck', 'cock','communication','sept','hiv', 'mouse', 'fever','sun','disorder','behaviour','immune', 'valve', 'shoulder', 'storm', 'mortality', 'angry','nightmare','muscle','vitamin','chicago','ear','injury','strength','wolf','alcohol','tiger', 'panel','shape','cbs','horn', 'channel', 'bull', 'mouth','amsterdam','snake', 'surgery','canadian','korean','cattle','dick','protective','thumb', 'coal','lemon','smoke','baltimore','coffee','spring', 'tooth','guide', 'fear', 'sheet','leather','danger','boston', 'illness','minnesota', 'bridge', 'arrest','virgin','egg','toxic', 'crowd', 'milk', 'asia','spin', 'anxiety','bloody','nervous','dental','wood','silent','circulation', 'barrier','migration', 'sand','mobility','divide', 'sweat','wave','harvest','corn','vessel','founder', 'tobacco','athletic', 'drain', 'surface', 'singapore','recovery', 'silk', 'tongue','tail', 'discrimination','react', 'pollution', 'tree', 'irish','map', 'polar', 'sick', 'maintenance','grow','guinea','stem','environmental', 'lung','organ', 'chicken', 'lesbian','transport','feedback', 'frequency','neck','spanish', 'therapy', 'hunter','mac','oak', 'tide','intellectual', 'moscow','portuguese', 'clinic','colorado','attraction','determination','iron','bend','drug', 'knee','battle', 'sugar','pocket', 'bear','russian','motor', 'liver','cigarette', 'musical', 'delivery','bat','gay', 'recognition','clinical', 'intelligent', 'bubble','testosterone','mutation', 'aggression','grape','arthritis','infectious','obese', 'maternity', 'hereditary','crab', 'landmark', 'shalt', 'plasma', 'endemic','sensory', 'zinc', 'inflammation', 'transitional','fatty','cone','juvenile', 'tolerant','sponge','veterinary','agony','insect', 'platinum','distress','colon','cardiac','homo','node','bleed', 'heel', 'climax', 'homicide','goat','rift', 'cerebral','claw', 'ape','vegan','lethal','feminine','schizophrenia','tor','enzyme', 'hormone','tomato', 'allergic', 'penetration','fork','nutrient', 'bipolar', 'progression','abdominal','heroin', 'gut','autopsy','masculine','bulb', 'differentiate','irrigation','headache','mast', 'rosemary', 'rib','oyster','suppression','dementia','paralysis', 'reproduction','sensor', 'cholesterol','flora', 'cinnamon','duct', 'maturity', 'muscular', 'autonomous','breadth','fisher','allergy','rag','fda','ark','porter', 'stimulus','maze','epidemic', 'bengal','lang', 'abdomen', 'hepatitis','displacement','carrot', 'ecology','abduction', 'heartbeat','parasite','algeria','deaf', 'psychiatric', 'cardiovascular', 'elbow', 'medicare','dwarf','vaccine','malaria','psa','titanium','substitution', 'fatigue', 'trout', 'maxwell','fiber','hoover', 'behavioral', 'tidal','metabolic','obesity', 'sensitivity','calcium','vein', 'aluminum','turtle','asteroid','poisonous', 'sperm','cox', 'adolescent','urine','contagious','gluten','ankle','butterfly', 'toxicity', 'diarrhea', 'crossover','fetus','floral','iris','hostility','collateral', 'confrontation', 'rehab', 'mentality', 'vaccination', 'clearance', 'salvage','anal','penetrate', 'prostate', 'repression', 'webster', 'protector', 'psychic', 'projection', 'inheritance','dentist', 'audition','disposition', 'vinegar', 'adaptive', 'ecstasy', 'appetite','instability', 'restless', 'swan','locomotive', 'erection', 'mushroom','antibiotic','owl','invasive', 'pneumonia', 'pyramid', 'aquatic', 'plantation', 'curtis', 'junction', 'neural', 'eruption', 'insanity','homosexual','caffeine','promoter','vomit','bloom', 'asylum', 'homosexuality', 'onion','hairy','flavor', 'regulator', 'membrane', 'immature', 'oxide', 'sting','bind', 'skeleton','amp','developmental', 'swell', 'plague','spinal', 'cbc','goose', 'flu', 'peripheral', 'motive','adapter', 'fog', 'alcoholic','spine','nap','reg','insulin', 'larvae','glucose', 'stationary', 'autism','tumor','mole','envy','tuberculosis','starvation', 'segregation','frog', 'immunity','nucleus','respiratory', 'resistant','inflammatory','fusion', 'venom', 'banana','mosquito','reproductive','proliferation','digestive','cohesion','deprivation', 'beetle','serpent', 'obstruction','insomnia', 'blindness','nasal', 'leukemia','amnesia','fungus','diabetic','snail','lac','moth', 'peacock','colonization','embryo', 'cortex', 'pediatric','vascular', 'magnesium', 'cholera', 'fetal','exhaustion','fungi', 'intolerance','marrow','innate', 'acne','mule', 'provocative','semen', 'cervical','nicotine','retina', 'cdc','thyroid','influenza', 'differentiation','alcoholism','vaginal','benign','adrenaline','pulmonary', 'migrate','atp','perennial']
    return keepterms

def generate_stop_grams():
    stoponegram = ['fulli','arise','luciferase','clearly','mechanistic','elevate','fxr','atypycal','orfs','spectroscopic','orthology','strikingly','reciprocal','devoid','transrciptomics','genomics','in-silico','p005', 'p00001', 'p0001','ass','subsequent','manner','aid','score','correct','gain','receive','partner','neither','temperature','physiological','induction','hybrid','crucial','similarity','directly','phase','core','impair','define','specificity','double','mass','kda','basis','marker','interestingly','regulatory','affinity','nuclear','recently','screen','fold','structural','appear','characterization','critical','library','direct', 'independent', 'approximately', 'electron','consist','band','primary','conclude','separate','stimulate','block','alpha','beta','larger','smaller','additional','chain','carry','secondary','pre','label','internal','hybridization','transfer', 'alter','accumulate', 'restriction','modification','frame','chromosome','extra','eukaryote', 'synthesis', 'rpl22', 'fragment', 'plasmid', 'precursor', 'bacterial','derive', 'contain', 'yeast', 'amino', 'peptide', 'transcribe', 'initiation','bp', 'residue', 'clone', 'helix','locate', 'infection', 'terminus', 'mutant', 'terminal', 'bacteria','virus', 'strand', 'pair', 'microscopy', 'gel', 'reaction', 'contrast', 'rnase', 'fraction', 'host', 'repeat', 'iii', 'generate', 'probe', 'tissue', 'electrophoresis', 'domain', 'obtain', 'formation', 'eukaryotic', 'purify', 'cleavage','pattern', 'characteristic', 'clade', 'identical', 'product','transcription','complex','proteomic','relevance', 'additionally', 'derivative', 'functionality', 'differential', 'clue', 'insight', 'involvement', 'stability', 'pcr', 'spp', 'favor', 'gc', 'conserve', 'cod', 'amplification','discrete', 'complementary', 'mix', 'correspond', 'transcript', 'abundant', 'divergent', 'comparative', 'vicinity', 'interact', 'circular', 'rna', 'interference', 'stag', 'remarkably', 'withdrawal', 'notably', 'measurement', 'mature', 'rely', 'distinctive', 'expand', 'newly', 'strain', 'participate', 'aspect', 'meanwhile', 'barely', 'discovery', 'formula', 'switch', 'bind', 'closely', 'genome', 'display', 'raise', 'toward', 'apply', 'copy', 'mostly', 'prevent', 'progress', 'animal', 'reference', 'indeed', 'location', 'responsible', 'relevant', 'throughout', 'specifically','element', 'v79','structure','absent', 'regulation', 'kb','additionally','transcriptomics','transcriptomics','ribonomics','ribonomic','constitutive','encode','relevance','seq','pcr','electrophoretic','inhibit','mediate','activate','mammalian','inhibition', 'regulate','inhibitor','activation','cell','agent','express','molecule','downstream','activation', 'inhibitor', 'regulate','inhibition', 'activate', 'assay','inhibit', 'treat', 'promote', 'anti', 'mediate', 'suppress', 'cellular',"addres","identification","poorly","confirm","extent","strongly","distinct","moreover","experimental","mrna","cdna","overexpression","inhibited","attenuated","transcriptional","knockdown","proccesing","intracellular","dataset","polymorphism","phenotype","subunit","mediate","understand","wether","constitutively","homologues","determinant","transiently","blot","recombinant","modulate","nucleotide","polypeptide","putative","homologous","exon","elucidate","concomitant","homolog","chromosomal","loci","biomarker","causative","overexpressed","recombination","inducible","mechanistically","ubiquitously","correlate","subtype","determinant","perturbation","implicate","cluster","probands","counteract","overlap","supplementation","alternation","underrepresented","transfected","upregulated","homology","differentially","uncharacterized","upregulation","quantitation","colocalization","localization","orthologue","ortholog","downregulated","category","composition","primer","linkage","conservation","diversity","canonical","administration","intrinsic","structurally","soluble","comprenhensive","availability","adaptation","confer","highest","reactive","pool","prepare","vector","situ","classical","classification","coordinate","independently","positively","null","localize","specie","termini","deduce","effector","discus","activity","require","thousand","protein","typically","occur","essential","demonstrate","component","active","initial","reduce","train","context",'μmol', 'μm', 'μg', 'µm', 'yr', 'yield', 'wo','widespread', 'widely', 'wide', 'whereas', 'weak', 'volume', 'vivo', 'vitro','via', 'versus', 'verify', 'vary', 'variety', 'variation', 'variant', 'variance', 'variable', 'variability', 'valuable','validate','utilize', 'utilization','utility', 'useful','uptake', 'upper', 'update', 'unlike', 'unknown', 'unite', 'unit', 'unique', 'underlie', 'undergo', 'unclear', 'uncertainty', 'un', 'ultimately', 'ubiquitous','sub', 'sup', 'gene', 'suggest', 'sample', 'analysis', 'identify', 'associate', 'indicate', 'variation', 'concentration', 'impact', 'genetic', 'relate', 'distribution', 'compare', 'affect', 'factor', 'sequence', 'sp', 'condition', 'scale', 'describe', 'reveal', 'method', 'interaction', 'significantly', 'investigate', 'influence', 'mechanism', 'acid', 'function', 'observe', 'genus', 'approach', 'significant', 'abundance', 'decrease', 'analyse','expression', 'determine', 'measure', 'limit', 'examine', 'relative', 'remain', 'collect', 'induce', 'experiment', 'strategy', 'molecular', 'estimate', 'sit', 'target', 'respectively', 'nutrient', 'dynamic', 'content', 'protect', 'evolution', 'wild', 'global', 'predict','multiple', 'presence','dna', 'resource', 'develop', 'organism', 'taxa', 'novel', 'positive', 'consider', 'evolutionary', 'contribute', 'improve', 'phylogenetic', 'signal', 'native', 'evaluate', 'addition', 'ratio', 'produce', 'hypothesis','involve', 'represent', 'characterize', 'enhance','highly', 'detect', 'focus', 'functional', 'importance', 'compound', 'overall', 'differ', 'le', 'greater', 'highlight', 'negative', 'copyright', 'combine', 'lineage','biological', 'exhibit','reduction', 'correlation', 'occurrence', 'consistent', 'quantify', 'aim', 'combination','conduct', 'analyze', 'pathway','lack', 'ability', 'explain','challenge', 'flow', 'effective', 'previously','perform', 'propose', 'isolate','explore', 'transmission','maintain', 'length', 'despite','establish', 'stable', 'facilitate', 'furthermore','dependent', 'enzyme', 'consequence', 'potentially', 'mg', 'feature', 'profile', 'quantitative','select', 'monitor', 'tool','relatively','comparison', 'link', 'framework', 'distance', 'evolve','genomic']
    stoptwogram = ['v v', 'control region','reverse transcriptase','base pairing',  '1 2','two different','x1 x2','double stranded','x2 x3', 'free system', 'stringent control', 'ph 8', 'cryo em', 'e dispar','mass spectrometry', '1 10', 'also observed', 'data show','relationships among', 'decision making','0 01', 'human health', 'pancreatic cancer','95 ci', 'ci 1', 'better understanding','within species', '0 3', 'well studied', 'also provide', 'species within', 'decision makers','rights reserved','p 0', 'land use', 'results showed', 'poorly understood', 'present study', 'two new', 'results show', 'life history', 'two species','0 001','time spent', 'results indicated', '0 5', 'results provide', 'trade offs', '5 ht', 'g c', '3 5', '1 5', '2 5', 'nov n', 'may also', '0 1', 'higher levels', '1 0', 'better understand', 'c n', 'study demonstrates', 'also found', 'times higher','well understood', 'bacterial communities', 'years ago', 'future studies', '1 1', 'scanning electron',"suggest","results","demonstrate","functionally","distinct","understood","data","play","conserved","members","also","study","reported","simultaneously","analysis","allowed","characterize",'important role',"gene expression","cell lines","expression levels","transcription factor","mrna expression","closely related","mammalian cells","family genes","target genes","amino acid","amino acids","differentially expressed","important roles","protein expression","taken together","significantly increased","first time","previous studies","acid sequence",'0 05', 'real time', 'well known', 'long term', 'large number', 'free energy', 'type 1', 'open reading', 'quality control', 'next generation','copy number', 'type ii', 'cell free', 'widely used', 'high resolution', 'provide evidence', 'e g', 'cell line', 'full length', 'large scale','cell type', 'wild type', 'recent studies', 'wide range', 'expression level', 'parameter','et al', 'steady state', 'entry site', 'three dimensional', 'cell types', 'reading frame', 'host cell', 'highly expressed', 'protein levels', 'sub 2', '2 sub', 'high affinity', 'genes involved','largely unknown', 'related genes', 'å resolution', 'base pairs','proteins involved', 'secondary structures', 'small molecule', 'recent advances', 'non coding', 'reading frames', 'single stranded', 'single molecule', 'expression patterns', 'dependent manner', 'biological processes', 'protein protein', 'protein 1', 'crystal structures', 'high throughput', 'binding sites', 'de novo', 'c virus', 'remains unclear', 'expressed genes', 'terminal domain', 'genome wide', 'cellular functions', 'cellular processes', 'green fluorescent', 'molecular mechanism', 'mechanisms underlying', 'protein coding', 'molecular mechanisms', 'coding sequence', 'protein folding' ,'molecular dynamics', 'transmission electron', 'r proteins','sup 2', '2 sup', 'sup sup', 'nucleic acids','fluorescent protein','genes encoding','coding rnas','western blotting', 'site directed','x ray',]
    stopthreegram = ['lc ms ms','','x1 x2 x3', 'χ 2 p', '2 p 0','coding genes 22','degrees c',"amino acid sequence",'p 0 01','p 0 05',"p 0 0001","p 0 001",'sup 2 sup','open reading frames','sub 2 sub','green fluorescent protein', 'c virus hcv']
    stopgram = [stoponegram,stoptwogram,stopthreegram]
    return stopgram


def cleanwords(worddict,stopgram,keepterms):
    stemmer = SnowballStemmer("english")
    lemmatizer = nltk.stem.WordNetLemmatizer()
    stopwords = nltk.corpus.stopwords.words('english') #stopword
    listarticles = list(worddict.keys())
    for z in stopgram[0]:
        stopwords.append(z) # agrego a los stopwords las relevantes a los onegram
    for article in listarticles:
        cleaned_words = [w for w in worddict[article] if (len(w) > 2 and not w[0].isdigit() and not w[-1] == "-")]
        cleaned_words = [w for w in cleaned_words if w not in stopwords]
        cleaned_words = [lemmatizer.lemmatize(l,pos="v") for l in cleaned_words]
        cleaned_words = [lemmatizer.lemmatize(l) if l[-1] == "s" else l for l in cleaned_words]
        cleaned_words = [stemmer.stem(l) if (l.endswith("ing") or l.endswith("lly")) else l for l in cleaned_words]
        cleaned_words = [w for w in cleaned_words if w not in stopwords and (zipf_frequency(w, 'en') <= 3.4 or w in keepterms)]
        cleaned_words = [w for w in cleaned_words if (len(w) > 2)]
        worddict[article] = cleaned_words
    return worddict


#Tiempos que demora:


def dict_gram(word_dict,mostcommon=1000):
    gramdict = defaultdict(dict)
    bag_of_words = list(itertools.chain.from_iterable(word_dict.values()))
    freq = nltk.FreqDist(bag_of_words)
    for x in freq.most_common(mostcommon):
        term = x[0]
        freq_value = x[1]
        id_abstracts = [i for i,y in word_dict.items() if term in set(y)]
        narticles = len(id_abstracts)
        if narticles > 0:
            gramdict[term] = [freq_value,narticles,",".join(id_abstracts)] #numero de articulos en que aparece mencionado el grama
        else:
            pass
    return gramdict

def n_containing(word, bloblist):
    return sum(1 for blob in bloblist if word in blob.words)

def idf(word, bloblist):
    return math.log(len(bloblist) / (1 + n_containing(word, bloblist)))

def idf_dict(worddict,dictgram):
    idfscores = {}
    bloblist = [tb(" ".join(x)) for x in worddict]
    for term in dictgram.keys():
        idfscores[term] = idf(term,bloblist)
    return idfscores

def get_gene_term_dicts(gramdict,idfscore):
    termdict = defaultdict(dict)
    geneterms = defaultdict(dict)
    geneNterms = defaultdict(dict)
    gramlist = gramdict.keys()
    for term in gramlist:
        for index,gene in enumerate(genepmids):
            geneterms[term][gene] = ",".join([x for x in gramdict[term][2].split(",") if x in genepmids[gene].split(",")])
            if geneterms[term][gene] != "":
                geneNterms[term][gene] = geneterms[term][gene].count(",") +1 
            else:
                geneNterms[term][gene] = 0
        termdict[term]["Articles"] = gramdict[term][1] #Numero de articulos que mencionan el termino 
        termdict[term]["Genes"] = sum([1 for gene in geneterms[term] if geneterms[term][gene] != ""]) #Numero de genes que mencionan el termino 
        termdict[term]["Freq_Term"] = int(gramdict[term][0]) #Numero de veces que se menciona el termino por articulo que lo hace. 
        termdict[term]["Zip_Score"] = round(zipf_frequency(term, 'en'),3) #zipf score (general freq of the word)
        termdict[term]["IDF_Score"] = round(idfscore[term],3)
        termdict[term]["Publications"] = gramdict[term][2]
    return  termdict,geneterms, geneNterms   

def create_co_occurrence_matrix(genepmids,gramdict,geneNterms,outputfile="matrix.csv"):
    dfco = pd.DataFrame(columns = list(genepmids.keys()), index = gramdict.keys())
    dfco[:] = int(0)
    for gene in list(genepmids.keys()):
        for term in gramdict:
            dfco[gene][term] = geneNterms[term][gene]
    dfco.to_csv(outputfile,index=True,header=True,sep=",")

def create_terms_info_dataframe(termdict,outputfile="term-info.csv"):
    dfterms = pd.DataFrame(termdict).T
    dfterms["Terms"] = dfterms.index.tolist()
    dfterms = dfterms[["Terms","Genes","Zip_Score","IDF_Score","Articles","Freq_Term","Publications"]]
    dfterms.to_csv(outputfile,index=False,header=True,sep=",")


def create_gene_artciles_dataframe(geneterms,outputfile="Publications-by_gene.csv"):
    dfgenesbyarticles =pd.DataFrame(geneterms).T
    dfgenesbyarticles["Terms"] = dfgenesbyarticles.index.tolist()
    dfgenesbyarticles.to_csv(outputfile,index=False,header=True,sep=",")


#Primer paso: Cargo los archivos resultantes de paperBLAST. Los que me indican los genes por articulo.

data_folder = Path("C:/Users/santi/") #Lo hago asi porque es windows....
file_to_open = data_folder / "Articles.tab" #Lo hago asi porque es windows....
df = pd.read_csv(file_to_open,sep="\t",header=0) #Cuento con 472 genes...
df.dropna(inplace=True) #Para uno de ellos ningun articulo paso los filtros
genepmids = ids_by_gene(df[:100]) #Gurado esa informacion en un diccionario.

#Segundo paso obtengo una lista de abstracts y sus ids correspondientes

start_time = time.time()
abstractdict = genes_list_abstarcts(genepmids) # 38 segundos para 100 genes obteniendo 8489
elapsed_time = time.time() - start_time

print(elapsed_time," <- obtener la lista de articulos: 5 genes")

#Obtengo una lista de palabras por articulo
start_time = time.time()
word_dict =get_words(abstractdict)
elapsed_time = time.time() - start_time

print(elapsed_time," <- obtener la lista de palabras por articulos: 5 genes")

# Guardo informacion del los stop words y de los keepwords
stopgram = generate_stop_grams() #palabras a eliminar.
keepterms = get_keepterms()

#Obtengo palabras mas limpias

start_time = time.time()
word_dict = cleanwords(word_dict,stopgram,keepterms)
elapsed_time = time.time() - start_time

print(elapsed_time," <- obtener la lista de palabras limpias por articulos: 5 genes")

#Obtengo el diccionario de terminos por articulo: 297 segundos

start_time = time.time()
dictgram = dict_gram(word_dict,mostcommon=10000)
elapsed_time = time.time() - start_time

print(elapsed_time," <- obtener diccionario de terminos (10000) por articulo: 5 genes")

#Demora 33 segundo en calcular el IDF para todos los terminos
start_time = time.time()
idfscore = idf_dict(word_dict,dictgram)
elapsed_time = time.time() - start_time

print(elapsed_time," <- Calculo del IDF score para los 10000 terminos: 5 genes")


#Obtengo los diccionarios para realizar las visualizaciones: 150 segundos

start_time = time.time()
termdict,geneterms, geneNterms = get_gene_term_dicts(dictgram,idfscore)
elapsed_time = time.time() - start_time

print(elapsed_time," <- Obtengo dict terminos y genes: 5 genes")


#Creo matriz de co-ocurrencia. Termino y gen. 76 segundos
start_time = time.time()
create_co_occurrence_matrix(genepmids,dictgram,geneNterms,outputfile="matrix.csv")
elapsed_time = time.time() - start_time

print(elapsed_time," <- Matriz de co-ocurrencia 10000 terminos: 5 genes")


#Creo dataframe de informacion de termino. 1 segundo
start_time = time.time()
create_terms_info_dataframe(termdict)
elapsed_time = time.time() - start_time

print(elapsed_time," <- Informacion de los 10000 terminos: 5 genes")


#Creo dataframe de articulos por gen.
start_time = time.time()
create_gene_artciles_dataframe(geneterms)

elapsed_time = time.time() - start_time

print(elapsed_time," <- Informacion de articulo por gen: 5 genes")


# def getgram(n,abstractslist,idlist,wordsbyabstracts,stopgram):
#     gramdict = defaultdict(list)
#     lemmatizer = nltk.stem.WordNetLemmatizer()
#     stopwords = nltk.corpus.stopwords.words('english') #stopwords
#     abstracts = abstractslist
#     for z in stopgram[0]:
#         stopwords.append(z) # agrego a los stopwords las relevantes a los onegram
#     if n == 1:
#         cleaned_words = [[w for w in subabstract if (len(w) > 2 and not w[0].isdigit())] for subabstract in wordsbyabstracts]
#         cleaned_words = [[w for w in subabstract if w not in stopwords] for subabstract in cleaned_words]
#         cleaned_words = [[lemmatizer.lemmatize(l,pos="v") for l in subabstract] for subabstract in cleaned_words]
#         cleaned_words = [[lemmatizer.lemmatize(l) if l[-1] == "s" else l for l in subabstract]for subabstract in cleaned_words]
#         cleaned_words = [[stemmer.stem(l) if (l.endswith("ing") or l.endswith("lly") and zipf_frequency(l,"en") <= 4 ) else l for l in subabstract]for subabstract in cleaned_words]
#         cleaned_words = [[w for w in subabstract if w not in stopwords and zipf_frequency(w, 'en') <= 5]for subabstract in cleaned_words]
#         bag_of_words = list(itertools.chain.from_iterable(cleaned_words))
#         freq = nltk.FreqDist(bag_of_words)
#         freq_common = [x for x in freq.items() if x[1] >= 2]
#         for x in freq_common:
#             term = x[0]
#             freq_value = x[1]
#             in_abstract = [1 if re.search(r'\b' + term + r'\b', " ".join(y)) else 0 for y in cleaned_words]
#             narticles = sum(in_abstract)
#             id_abstracts = ",".join([v for i,v in enumerate(idlist) if  in_abstract[i] == 1])
#             if narticles > 0:
#                 gramdict[term] = [freq_value,narticles,id_abstracts] #numero de articulos en que aparece mencionado el grama
#             else:
#                 pass
#     elif n == 2:
#         bigrams = [[" ".join(map(str,b)) for b in nltk.bigrams(subabstract) if (b[0] not in stopwords and b[1] not in stopwords and " ".join(map(str,b)) not in stopgram[1])] for subabstract in wordsbyabstracts]
#         bag_of_words = list(itertools.chain.from_iterable(bigrams))
#         freq = nltk.FreqDist(bag_of_words)
#         freq_common = [x for x in freq.items() if x[1] >= len(abstracts)*0.01]
#         for x in freq_common:
#             term = x[0]
#             freq_value = x[1]
#             in_abstract = [1 if re.search(r'\b' + term + r'\b', " - ".join(y)) else 0 for y in bigrams]
#             narticles = sum(in_abstract)
#             id_abstracts = ",".join([v for i,v in enumerate(idlist) if  in_abstract[i] == 1])
#             if narticles > 0:
#                 gramdict[term] = [freq_value,narticles,id_abstracts] #numero de articulos en que aparece mencionado el grama
#             else:
#                 pass
#     elif n == 3:
#         trigrams = [[" ".join(map(str,b)) for b in nltk.trigrams(subabstract) if (b[0] not in stopwords and b[1] not in stopwords and b[2] not in stopwords and " ".join(map(str,b)) not in stopgram[2])] for subabstract in wordsbyabstracts]
#         bag_of_words = list(itertools.chain.from_iterable(trigrams))
#         freq = nltk.FreqDist(bag_of_words)
#         freq_common = [x for x in freq.items() if x[1] >= len(abstracts)*0.01]
#         for x in freq_common:
#             term = x[0]
#             freq_value = x[1]
#             in_abstract = [1 if re.search(r'\b' + term + r'\b', " - ".join(y)) else 0 for y in trigrams]
#             narticles = sum(in_abstract)
#             id_abstracts = ",".join([v for i,v in enumerate(idlist) if  in_abstract[i] == 1])
#             if narticles > 0:
#                 gramdict[term] = [freq_value,narticles,id_abstracts] #numero de articulos en que aparece mencionado el grama
#             else:
#                 pass
#     return gramdict

# def genes_info(genepmids,stopgram):
#     """
#     Function to obtain the one and bigrams associated with a specific gene. 
#     input: pmids: list of pubmed id. gene: gene who is associated with the pmids.
#     output: genedict: dictionary of genes and his associated one and bigrams.
#     """
#     full_abstracts = []
#     full_pubmedid = []
#     for gene,pmids in genepmids.items():
#         abstarct,pubmedid = get_abstarcts(pmids,gene)
#         full_abstracts += abstarct
#         full_pubmedid += pubmedid
#     wordsbyabstracts = get_words(full_abstracts)
#     onegramdict = getgram(1,full_abstracts,full_pubmedid,wordsbyabstracts,stopgram)
#     #twogramdict = getgram(2,full_abstracts,full_pubmedid,wordsbyabstracts,stopgram)
#     #threegramdict = getgram(3,full_abstracts,full_pubmedid,wordsbyabstracts,stopgram)
#     return onegramdict,full_abstracts,full_pubmedid

# def term_counts(genepmids,gramdict,full_abstracts,stopgram,type):
#     """
#     term_counts define the number of genes and the number of articles that a gram is present. Also count the number of genes which has articles that have the term.
#     input: gramdict: dictionary of terms and count of term (number of article where it appear) for a specific gene.
#     output: dictionary of terms that indicate the total gene which the term has connection with it 
#     """
#     zipterm = ""
#     termdict = defaultdict(list)
#     gramlist = [x for x in gramdict.keys()]
#     if type == 1:
#         idfscore = idf_dict(gramdict,gramlist,full_abstracts,stopgram)
#     for term in gramlist:
#         termdict[term] = [0]*(len(genepmids.keys())+5)
#         for index,gene in enumerate(genepmids):
#             narticlesbygene = sum([1 if x in  genepmids[gene].split(",") else 0 for x in gramdict[term][2].split(",")])
#             termdict[term][index+1] = narticlesbygene
#         termdict[term][0] = len(termdict[term]) - termdict[term].count(0) #number of genes associated with the term
#         termdict[term][-1] = gramdict[term][1] # number of articles that contain the term
#         termdict[term][-2] = gramdict[term][0] # frequency of the term in the document
#         zipterm = term.split(" ")
#         flag = any(1 if len(x) >= 2 else 0 for x in zipterm)
#         if flag:
#             zipterm = [x for x in zipterm if len(x) >= 2 ]
#         termdict[term][-3] = round(sum(zipf_frequency(t, 'en') for t in zipterm)/len(zipterm),3) #zipf score (general freq of the word)
#         if type == 1:
#             termdict[term][-4] = round(idfscore[term],3) #idf score
#         else:
#             termdict[term][-4] = 0 #idf score only computed for one grams.
#     return termdict

# def compute_score(termdict,genepmids):
#     sorteddict = sorted(termdict.items(), key=lambda t: t[1][0],reverse=True)
#     termscore = {}
#     for k,v in sorteddict:
#         plus = [1 if len(k) < 4 else 0]
#         highfreq = [1 if (v[-2]/v[-1]) > 0 else 0]
#         score = ((v[0]/len(genepmids))*10 - 2*max(v[1:-4])/v[-1] -1.5*(v[-3]) -(len(genepmids)/v[-2]) + plus[0] -v[-4] + highfreq[0]) + 20
#         termscore[k] = round(score,3)
#     return termscore

# def abstractwordcloud(termscore,pic,name):
#     print("Creating WordCloud!")
#     maskpic = np.array(Image.open(pic))
#     image_colors = ImageColorGenerator(maskpic)
#     wc = WordCloud(width=900,height=500, max_words=len(termscore),relative_scaling=1,normalize_plurals=False, mask=maskpic,random_state=1,margin=5).generate_from_frequencies(termscore)
#     wc.to_file(name + "_WordCloud.png")
#     plt.imshow(wc,interpolation="bilinear")
#     plt.axis("off")
#     plt.figure()
#     plt.close()
#     return None

# # Funciones Auxiliares

# def freqterms(gramdict,nterminf,ntermsup):
#     sortedict = sorted(gramdict.items(), key=lambda t: t[1][0],reverse=True)
#     count = 0
#     terms = []
#     for k,v in sortedict:
#         count +=1
#         if count > nterminf and count < ntermsup:
#             terms.append(k)
#     return terms

# # ##### MAIN FUNCTIONS #####

# data_folder = Path("C:/Users/santi/")

# file_to_open = data_folder / "Articles.tab"

# df = pd.read_csv(file_to_open,sep="\t",header=0)

# ############### 

# stopgram = generate_stop_grams()
# df.dropna(inplace=True)
# genepmids = ids_by_gene(df)
# onegramdict,full_abstracts,full_pubmedid = genes_info(genepmids,stopgram)
# #termdict = term_counts(genepmids,onegramdict,full_abstracts,stopgram,1)
# #termscore = compute_score(termdict,genepmids)
# #twotermdict = term_counts(genepmids,twogramdict,full_abstracts,stopgram,2) 
# #abstractwordcloud(termscore,"brain-990x622.png","guille")









