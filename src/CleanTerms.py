from src.HelperGrams import generate_stop_grams, get_keepterms, nltk
from nltk.stem.snowball import SnowballStemmer
from nltk import FreqDist
import itertools
from wordfreq import zipf_frequency
from collections import defaultdict
import string
import logging


# Download nltk data

nltk.download("wordnet")

nltk.download("punkt")

nltk.download("stopwords")

def get_words(abstractdict):
    """ Tokenization of articles abstracts. Breaking up sequences of strings into pieces (words/terms). For each abstract remove punctuation and transform to lower case. 
    
    Arguments:
        abstractdict {[dict]} -- Dictionary of articles, where the key is the pubmed id and the value is the abstract 
    
    Returns:
        [type] -- Dictionary of articles, where the key is the pubmed id and the value is a list of tokens.  
    """

    worddict = {} #empty dict
    removepunct = string.punctuation.replace("-","").replace("_","") #keep - and _ in text other punctuation removed
    table = str.maketrans('', '', removepunct) # create table to make the replacements
    for article in abstractdict: 
        worddict[article] = ([*map(str.lower, nltk.word_tokenize(abstractdict[article].translate(table)))]) #Tokenization and pre-processing
    return worddict

def cleanwords(abstractdict,keep,remove,zipf):
    """ Transform the list of word of each articles into a cleaned list (words without common terms, and trying to reduce the amount of words via lemmatization and stemmatization)
    
    Arguments:
        abstractdict {[dict]} -- Dictionary of articles, where the key is the pubmed id and the value is the abstract
    
    Returns:
        [dict] --Dictionary of articles, where the key is the pubmed id and the value is a list of cleaned (remove common words, stemm, lemma, etc.) tokens 
    """
    stopgrams = generate_stop_grams()
    stopgrams += remove
    keepgrams = get_keepterms()
    keepgrams += keep
    worddict = get_words(abstractdict)
    stemmer = SnowballStemmer("english") #stemmer  
    lemmatizer = nltk.stem.WordNetLemmatizer() #lemmatizer
    for article in worddict.keys():
        cleaned_words = [w for w in worddict[article] if (len(w) > 2 and not w[0].isdigit() and not w[-1] == "-")] # remove tokens if len < 3, start with a number or end with -.
        cleaned_words = [w for w in cleaned_words if w not in stopgrams] # remove common words 
        cleaned_words = [lemmatizer.lemmatize(l,pos="v") for l in cleaned_words] #transform to verb form
        cleaned_words = [lemmatizer.lemmatize(l) if l[-1] == "s" else l for l in cleaned_words] #try to remove plurals from not s ending words.
        #cleaned_words = [stemmer.stem(l) if (l.endswith("ing") or l.endswith("lly")) else l for l in cleaned_words] #Dealing with words ending in ing and lly because lemmatizer don work too good.
        cleaned_words = [w for w in cleaned_words if (zipf_frequency(w, 'en') <= zipf or w in keepgrams)] #keep word if zipf score is low (not frequent in english) or if it is in keep terms (words common in biology))
        cleaned_words = [w for w in cleaned_words if (len(w) > 2 and not w in stopgrams)]
        worddict[article] = cleaned_words
    return worddict

def dict_gram(worddict,mostcommon=1000,minfreq=2):
    """ Create a dictionary of terms (key).Each term need to appear at least minfreq in the corpus.
    
    Arguments:
        worddict {[dict]} -- Dictionary of articles, where the key is the pubmed id and the value is a list of tokens
    
    Keyword Arguments:
        mostcommon {int} -- Number of most frequent terms to include in the analysis  (default: {1000})
        minfreq {int} -- Minimum freq needed to be considerer (default: {2})
    
    Returns:
        [dict] -- Dictionary of terms (key). Values are [freq,# articles, pubmed ids]
    """
    gramdict = defaultdict(dict) #create empty dictionary
    bag_of_words = list(itertools.chain.from_iterable(worddict.values())) #create a list of all words from all articles and it associated frequency.
    freq = FreqDist(bag_of_words)
    logging.info("Total terms %i" %(len(freq)))
    for x in freq.most_common(mostcommon): #keep only most common words
        term = x[0]
        freq_value = x[1]
        if freq_value >= minfreq:
            id_abstracts = [i for i,y in worddict.items() if term in set(y)]
            narticles = len(id_abstracts)
        else:
            continue
        if narticles > 0:
            gramdict[term] = [freq_value,narticles,",".join(id_abstracts)] #numero de articulos en que aparece mencionado el grama
        else:
            pass
    logging.info("Total terms to be analyzed %i" %(len(gramdict)))
    return gramdict

