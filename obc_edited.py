"""
"""

import xml.etree.ElementTree as ET

import requests
import datetime as dt
import os

OBC_IDE_URL = "http://devapi.onbc.io/publications"
path1=os.path.dirname(os.path.abspath(__file__))

class PubmedArticle():

    def __init__(self, xml):  # start from xml as string for now
        self.xml = xml
        data = self._parse_xml(xml)
        self.doi = data["doi"]
        self.pmid = data["pmid"]

    def _parse_xml(self, root):
        # parse the data we're interested from the xml (passed as a string)
        # return the data as a dictionary
        # function assumes that xml root is a PubMedArticle node,
        # not a PubmedArticleSet node
        result = {"pmid": None, "doi": None}
        for pubid in root.findall('PubmedArticle'):
	    print(pubid.findall('MedlineCitation/PMID'))
            for med in pubid.findall('MedlineCitation/PMID'):
                result["pmid"] = med.text
                print(med.text)
            for doi in pubid.findall('MedlineCitation/Article/ELocationID'):
                if med.get('EIdType')=='doi':
                    result["doi"] = doi.text
                   
        #for article_id in root.find("PubmedArticle ").findall("MedlineCitation"):
        #    if article_id.get("PMID") == "PMID":
        #        result["pmid"] = article_id.text
        #    elif article_id.get("EIdType") == "doi":
        #        result["doi"] = article_id.text
        #    else:
        #       continue

        return result

        # these below lines are possible, but more more complicated than just using <ArticleIdList>

        # if root.find("MedlineCitation").find("PMID"):
        #     result["pmid"] = root.find("MedlineCitation")[0].text
        # else:
        #     result["pmid"] = None

        # result["doi"] = None
        # for eid in root.find("MedlineCitation").find("Article").findall("ELocationID"):
        #     if eid.get("EIdType") == "doi":
        #         result["doi"] = eid.text

        #return result

    def publication_in_obc(self):
        """Return True if an article is found in the obc.ide /publications api response.

        Otherwise, returns False. Checks against pmid first, then doi if pmid is
        not found. Accepts one argument, a PubmedArticle instance.
        """
        entries = requests.get(OBC_IDE_URL).json()
        
        for entry in entries:
            #print(self.pmid)
            print(entry.get("pmid"))
            #print(self.doi)
            #print(entry.get("doi", None))
            if self.pmid == entry.get("pmid", None):

                #print("matched",self.pmid)
                return True
    
            else:
                if self.doi == entry.get("doi", None):
                    #print("mactched",self.doi)
                    return True
        return False


def sort_xml(path):
    # sort PubmedArticles from specified .xml path into two buckets:
    # those that are found in obc, and those that are not found in obc
    in_obc = []
    not_in_obc = []
    from collections import defaultdict
    paper = defaultdict(list)
    grant=[]
    # loop thru the <PubmedArticleSet> node
    tree = ET.parse(path)
    root = tree.getroot()
    for node in root:
        article = PubmedArticle(node)
        
        if article.publication_in_obc():
            #print("pMID",node.find('./MedlineCitation/PMID').text)
            #print(ET.tostring(node))
            in_obc.append(article)
            
            
        else:
            for l in (node.findall('./MedlineCitation/Article/GrantList/Grant/GrantID')):
                grantid= l.text.split(" ")
                #print(pubid)
                if len(grantid)==3:
                    grantid=grantid[1]+grantid[2]
                    paper["".join(grantid)].append(node.find('./MedlineCitation')[0].text)
                else:
                    paper[" ".join(grantid)].append(node.find('./MedlineCitation')[0].text)
            not_in_obc.append(article)
    
    return in_obc, not_in_obc,paper


def write_to_xml(lst, path):
    # given a list of [PubmedArticle], write the list to xml file
    strs = [ET.tostring(a.xml).decode('utf-8') + "\n" for a in lst]
    with open(path, 'w') as f:
        #f.writelines('<?xml version="1.0" ?>'+'\n')
        #f.writelines('<!DOCTYPE PubmedArticleSet PUBLIC "-//NLM//DTD PubMedArticle, 1st June 2018//EN" "https://dtd.nlm.nih.gov/ncbi/pubmed/out/pubmed_180601.dtd">'+'\n')
        f.writelines("<PubmedArticleSet>"+'\n')
        f.writelines(strs)
        f.writelines('\n'+"</PubmedArticleSet>")
    return


#def get_xml():
#    today=dt.datetime.today().strftime(("%Y-%m-%d"))
#    filepath=path1+"\\"+("most_recent_publications_{}.xml".format(today))
#    in_obc, not_in_obc,paper = sort_xml(filepath)
#    print("In OBC: " + str(len(in_obc)))
#    print("Not In OBC: " + str(len(not_in_obc)))
#    write_to_xml(in_obc, path1+"\\"+'in_obc.xml')
#    write_to_xml(not_in_obc, path1+"\\"+'not_in_obc.xml')
#    return paper
if __name__ == '__main__':
     today=dt.datetime.today().strftime(("%Y-%m-%d"))
     filepath='C:\\Users\\AKSHATA\\Desktop\\OPS_Job\\obc-pubmed-master\\most_recent_publications_2019-03-19.xml'
     in_obc, not_in_obc,paper = sort_xml(filepath)
     print("In OBC: " + str(len(in_obc)))
     print("Not In OBC: " + str(len(not_in_obc)))
     write_to_xml(in_obc, path1+"\\"+'in_obc.xml')
     write_to_xml(not_in_obc, path1+"\\"+'not_in_obc.xml')
#get_xml()
