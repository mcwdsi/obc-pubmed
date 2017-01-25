""" Written in Python 3.5

The purpose of this program is to search PubMed for publications that were funded by the MIDAS grant with serial
number 'GM070708'."""

from Bio import Entrez
import xml.etree.ElementTree as ET


# Prompts user for tab-delimited text file containing grant numbers
def grants_prompt():
    try:
        prompt = input("Please enter the name of a tab-delimited text file containing the grant numbers. \n"
                    "\n"
                    "Also, please be aware that the grant numbers must end with '[gr]' (without single quotes) in order for the program to process\n"
                    "them. \n"
                    "\n"
                    "Input: ")
        fhand = open(prompt)
        return fhand
    except IOError:
        print("Please ensure that you spelled the file name correctly and that you're in the right directory.")


# Search PubMed by Grant No. and return primary ID's
def search(query):
    Entrez.email = 'diller17@ufl.edu'
    handle = Entrez.esearch(db='pubmed',
                            term=query,
                            retmax=1000,
                            )
    results = Entrez.read(handle)
    return results


# Search PubMed by primary ID's to obtain DocSums for each publication
def fetch_doc_sum(id_list):
    ids = ','.join(id_list)
    Entrez.email = 'diller17@ufl.edu'
    handle = Entrez.efetch(db='pubmed',
                           id=ids,
                           retmode='xml',
                           version='2.0')
    records = handle.read()
    return records


# Looks for path to XML tag and pulls out the text therein
def checkXML(XML, path):
    if XML.find(path) is not None:
        if XML.find(path).text is not None:
            return XML.find(path).text
        else:
            return ""
    else:
        return ""


# Parses out authors' last names and initials of first/middle names
def parse_authors(authors):
    authorList = []
    for author in authors:
        authorList.append(checkXML(author, './LastName') + " "
                          + checkXML(author, './Initials'))
    return authorList


# Parses out grant numbers
def parse_grants(grants):
    grant_list = []
    for grant in grants:
        grant_list.append(checkXML(grant, './GrantID'))
    return grant_list


# Extracts article title, author list, journal, date of publication, DOI, and PMID.
def extract_info(PubMedXML):
    tree = ET.fromstring(PubMedXML)
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


# Executes functions
if __name__ == '__main__':
    raw = grants_prompt()

    all_ids = set()

    for line in raw:
        line=line.rstrip()
        results = search(line)
        id_list = results['IdList']

        for id in id_list:
            all_ids.add(id)

    raw.close()

    papers = fetch_doc_sum(all_ids)
    #info = extract_info(papers)

    with open("output.xml", "a") as xml:
        xml.write(papers)

    print("Return count: ", len(all_ids), " papers")
