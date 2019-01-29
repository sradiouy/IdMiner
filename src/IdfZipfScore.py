import math
from wordfreq import zipf_frequency
from textblob import TextBlob as tb


def n_containing(word, bloblist):
    """returns the number of documents containing word. A generator expression is passed to the sum() function.
    
    Arguments:
        word {[type]} -- [description]
        bloblist {[type]} -- [description]
    
    Returns:
        [integer] -- number of documents containing word
    """

    return sum(1 for blob in bloblist if word in blob.words)

def idf(word, bloblist):
    """ Computes "inverse document frequency" which measures how common a word is among all documents in bloblist. The more common a word is, the lower its idf. 
    
    Arguments:
        word {[string]} -- A term to be evaluated
        bloblist {[list]} -- A bloblist
    
    Returns:
        [int] -- an inverse document frequency factor is incorporated which diminishes the weight of terms that occur very frequently in the document set and increases the weight of terms that occur rarely.
    """

    return round(math.log(len(bloblist) / (1 + n_containing(word, bloblist))),3)

def idf_dict(worddict,dictgram):
    """[summary]
    
    Arguments:
        worddict {[dict]} -- Dictionary of articles, where the key is the pubmed id and the value is a list of cleaned tokens
        dictgram {[dict]} -- Dictionary of terms (key). Values are [freq,# articles, pubmed ids]
    
    Returns:
        [dict] -- Dictionary of idf score by term.
    """

    idfscores = {}
    bloblist = [tb(" ".join(x)) for x in worddict.values()]
    for term in dictgram.keys():
        idfscores[term] = idf(term,bloblist)
    return idfscores


def zipf_score(dictgram):
    """ Calculate the Zipf score. Zipf's law states that given some corpus of natural language utterances, the frequency of any word is inversely proportional to its rank in the frequency table. Thus the most frequent word will occur approximately twice as often as the second most frequent word, three times as often as the third most frequent word, etc.: the rank-frequency distribution is an inverse relation.
    
    Arguments:
        dictgram {[dict]} -- Dictionary of terms (key). Values are [freq,# articles, pubmed ids]
    
    Returns:
        [dict] -- Zipf score
    """

    zipfdict = {}
    for term in dictgram:
        zipfdict[term] = round(zipf_frequency(term, 'en'),3)
    return zipfdict

