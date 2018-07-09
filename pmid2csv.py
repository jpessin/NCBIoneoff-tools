#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Read file of PMIDs (whitespace seperated), Query PubMed, via NCBI eutils (esearch),
Write parsed results to CSV."""
# Depends
# python 3.5+
# Biopython
# bs4, BeautifulSoup
# lxml
import csv
import os
import sys
import logging
import argparse

from datetime import datetime
from bs4 import BeautifulSoup
from Bio import Entrez as ez

####################
# CLI supplied -- going to be dynamic while working on so using names
# on complection switch to argparse or other cliparse &or configs

cliparser = argparse.ArgumentParser(
    description="Get Affiliations from NCBI-PubMed as .csv"
)

cliparser.add_argument("-i", "--pmids", help="REQUIRED   File with PMIDS to parse")
cliparser.add_argument("-o", "--outfile", help="REQUIRED   output file name")
cliparser.add_argument("--apikey", default=None, help="OPTIONAL  users NCBI api key")
cliparser.add_argument(
    "--keyfile", default=None, help="OPTIONAL   file containing only user NCBI api-key (UTF-8)"
)
cliparser.add_argument("--logging", help="OPTIONAL   turn logging on and use file path (debugging)")

cliargs = cliparser.parse_args()
###################
STARTTIME = datetime.now().isoformat()
COMMANDLINE = sys.argv

INFILE = cliargs.pmids
OUTFILE = cliargs.outfile
DOLOG = cliargs.logging
KEYFILE = cliargs.keyfile
KEYARG = cliargs.apikey


def checkreadfile(readfile=INFILE):
    """"Check readfile exists."""
    if not os.path.isfile(readfile):
        raise FileNotFoundError("can not find: " + readfile)

def checkwritefile(outfile=OUTFILE):
    """Check output dir exists, and filename availible."""
    outfile = os.path.realpath(outfile)
    outdir = os.path.dirname(outfile)
    if os.path.isfile(outfile):
        raise FileExistsError(outfile + " aready exists")
    elif not os.path.isdir(outdir):
        raise NotADirectoryError(outdir + " Not found")
    else:
        return outdir, outfile

##################
# do path checking
##################

checkreadfile(INFILE)
OUTDIR, OUTFILE = checkwritefile(OUTFILE)

##################
# logging setup - debugging, comment out for user runs
###############

if DOLOG:
    LOGGING_FILE = os.path.join(DOLOG, 'debug.log')
    logging.basicConfig(filename=LOGGING_FILE, level=logging.INFO)
    logging.info(STARTTIME)


###################
# Entrez setup bits
###################
ez.email = "jacob.pessin@einstein.yu.edu" # tool builders email
# ez.tool = toolname # default biopython
if KEYARG:
    ez.api_key = KEYARG
elif KEYFILE:
    with open(KEYFILE, 'rt') as kfile:
        FILEKEY = kfile.read()
    ez.api_key = FILEKEY

######################
# Read the file -> set(pubmedids)
#####################


def linecleaner(readline: str) -> set:
    """Split lines and elements into elementset."""
     # is this redundent when used inside a with like below? (what about files with mixed types)
    spltline = readline.splitlines()
    elems = set()
    for line in spltline:
        elems.update(line.split())
    return elems


def elementcleaner(elem: str) -> set:
    """Clean individual PMID elements.
    Helper func for linereader."""
    elem = elem.casefold().lstrip("pmid:")
    try:
        int(elem)
    except ValueError:
        raise TypeError("expected interger got: " + elem + "\n check your inputs")
    elemlen = len(elem)
    if  elemlen >= 1 and elemlen <= 8:
        return elem
    else:
        raise ValueError(
            "PMID's are intergers of span 1 to 8, got: " + elem + " with " + len(elem)
        )


def linereader(infile: str) -> set:
    """Read a file of cleaned PMID's."""
    with open(infile, 'rt') as ifile:
        elems = set()
        for line in ifile:
            elems.update(
                linecleaner(line)
            )
    elems = {elementcleaner(el) for el in elems}
    return elems


def pubmed_efetch(idstring) -> str:
    """Wrapper function for efetch with pre-sets."""
    records = ez.efetch(
        db="pubmed",
        id=idstring,
        post=True,
        retstart="0",
        rettype="MEDLINE",
        retmode="xml",
        retmax="9999"
    )
    return records


################################ 12
# xml parsing with BeautifulSoup
################################

# Helper func's for BeautifulSoup
BAILSTRING = "Tortuse didn't teach us so well - BeautifulSoup issues - blame the coder"


def children_astup(local_record):
    """ take a beautifulsoup record-branch, and returns newline removed
        tuple(branch.child)"""
    tmp = (child for child in local_record.children if child != '\n')
    return tuple(tmp)


def make_iterablerecords(entrez_data, parsemethod='lxml-xml'):
    """Turn a raw PubmedArticle or PubmedArticleSet to tuple of PubmedArticles.

    sets/subsets may be under a parent element '[document]'
    """
    ncbisettypes = ['RecordSet', 'PubmedArticleSet']
    ncbidatatypes = ['DocumentSummary', 'PubmedArticle']
    core = BeautifulSoup(entrez_data, parsemethod)

    if 'document' in core.name:
        child0 = children_astup(core)
        for kid in child0:
            if kid.name in ncbisettypes:
                data_lst = children_astup(kid)
        if data_lst not in locals(): # possible TODO: could be elif in preciding
            for kid in child0:
                if kid.name in ncbidatatypes:
                    data_lst = tuple(kid)
        elif core.name in ncbidatatypes:
            data_lst = children_astup(core)
        else:
            data_lst = core  # possible FIXME - other branches return tuple - ?tuple-ize?
    else:
        print(BAILSTRING, file=sys.stderr)
        #  ~ sys.exit()
        data_lst = core
    return data_lst


def articletitle(article):
    """Get Article Title as String"""
    if article.ArticleTitle:
        title = article.ArticleTitle.string
    else:
        title = "ArticleTitle not found"
    return title


def authorname(author):
    """Get Author Name."""
    if author.LastName and author.LastName.string:
        last = author.LastName.string
    else:
        last = str()

    if author.ForeName and author.ForeName.string:
        name = ", ".join((last, author.ForeName.string))
    elif author.Initials and author.Initials.string:
        name = ", ".join((last, author.Initials.string))
    elif author.CollectiveName and author.CollectiveName.string:
        name = author.CollectiveName.string
    else:
        name = last

    return name


def authoraffil(author):
    """Get Authors Affiliations as list of Strings."""
    affils = author.find_all('Affiliation')
    affiliations = [aff.string for aff in affils]
    return affiliations


def recs2array(recs):
    """Parse Records Tuple(PubmedArticleXML) array.
    [
    [PMID, title, AuthorName, affil1, afill2]
    ]"""
    info = list()
    cnt = 0
    for art in recs:
        cnt = cnt + 1
        if art.PMID:
            pmid = art.PMID.string
        elif art.ArticleIdlistList:
            for atr in children_astup(RECS[2].ArticleIdList):
                if "pubmed" in atr.attrs.values():
                    pmid = atr.string
        else:
            pmid = "PMIDNotFound-article" + cnt

        title = articletitle(art)

        authorlist = art.AuthorList.find_all('Author')
        for author in authorlist:
            name = authorname(author)
            row = authoraffil(author)
            row.insert(0, name)
            row.insert(0, title)
            row.insert(0, pmid)

            info.append(row)

    return info


def csvdump(array, outfile):
    """Write Array to csv file."""
    with open(outfile, 'wt', encoding='utf-8', newline='') as ofile:
        writer = csv.writer(ofile, delimiter=',', quoting=csv.QUOTE_ALL)
        for row in array:
            writer.writerow(row)


#####################
### doing the thing
#####################
FILEELEMENTS = linereader(INFILE)
IDSTRINGS = ",".join(FILEELEMENTS)
RECORDS = pubmed_efetch(IDSTRINGS)
RECS = make_iterablerecords(RECORDS)
RECSARRAY = recs2array(RECS)
csvdump(RECSARRAY, OUTFILE)
