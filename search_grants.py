"""
Written in Python 3.5

The purpose of this script is to search PubMed for publications that were funded by MIDAS grants."""

import requests
import xml.etree.ElementTree as ET
import datetime as dt
from Bio import Entrez


def grants_prompt():
    """Prompts user for text file containing grant numbers."""

    try:
        prompt = input("Please enter the name of a text file containing the grant numbers. \n\n"
                    "Also, please be aware that the grant numbers must end with '[gr]' (without single quotes) in order for the program to process\n"
                    "them. \n\n"
                    "Input: ")
        fhand = open(prompt)
        return fhand
    except IOError:
        print("Please ensure that you spelled the file name correctly and that you're in the right directory.")


def search(query):
    """Search PubMed by grant number and return primary IDs."""

    Entrez.email = 'diller17@ufl.edu'
    handle = Entrez.esearch(db='pubmed',
                            term=query,
                            retmax=1000,
                            )
    results = Entrez.read(handle)
    return results


def fetch_doc_sum(id_list):
    """Search PubMed by primary IDs to obtain DocSums for each publication."""

    ids = ','.join(id_list)
    Entrez.email = 'diller17@ufl.edu'
    handle = Entrez.efetch(db='pubmed',
                           id=ids,
                           retmode='xml',
                           version='2.0')
    records = handle.read()
    return records


def checkXML(XML, path):
    """Looks for path to XML tag and pulls out the text therein."""

    if XML.find(path) is not None:
        if XML.find(path).text is not None:
            return XML.find(path).text
        else:
            return ""
    else:
        return ""


def parse_authors(authors):
    """Parses out authors' last names and initials of first/middle names."""

    authorList = list()

    for author in authors:
        authorList.append(checkXML(author, './LastName') + " "
                          + checkXML(author, './Initials'))
    return authorList


def parse_grants(grants):
    """Parses grant numbers from DocSum XML."""

    grant_list = list()

    for grant in grants:
        grant_list.append(checkXML(grant, './GrantID'))
    return grant_list


def extract_info(pub_med_XML):
    """Not currently used.

    Extracts article title, author list, journal, date of publication, DOI, and PMID from the DocSum XML."""

    tree = ET.fromstring(pub_med_XML)
    articles = tree.findall('./PubmedArticle/MedlineCitation')

    papers = list()

    for article in articles:
        paper = dict()
        paper['TITLE'] = article.find('./Article/ArticleTitle').text
        paper['AUTHORS'] = parse_authors(article.findall('./Article/AuthorList/Author'))
        paper['PMID'] = article.find('./PMID').text
        paper['DOI'] = checkXML(article, './Article/ELocationID')
        paper['GRANT'] = parse_grants(article.findall('./Article/GrantList/Grant'))
        paper['JOURNAL'] = article.find('./Article/Journal/ISOAbbreviation').text

        papers.append(paper)

    return papers


obc_api_url = "http://api.onbc.io/publications"

# Executes functions
if __name__ == '__main__':
    grants_in = grants_prompt()

    all_ids = set()

    for line in grants_in:
        line = line.rstrip()
        results = search(line)
        id_list = results['IdList']

        for id in id_list:
            all_ids.add(id)

    grants_in.close()

    pmid_lst = list()

    in_obc = requests.get(obc_api_url).json()

    for publication in in_obc:
        pmid_lst.append(publication.get("pmid"))

    new_pub_ids = []

    for id in all_ids:
        if id not in pmid_lst:
            new_pub_ids.append(id)
        else : continue

    papers = fetch_doc_sum(new_pub_ids)
    #info = extract_info(papers)

    today = dt.datetime.today().strftime("%Y-%m-%d")
    fname = "most_recent_publications_{}.xml".format(today)

    with open(fname, "a") as xml_out:
        xml_out.write(papers)

    print("Number of new publications: ", len(new_pub_ids))

    #for id in new_pub_ids:
    #    print(id)
