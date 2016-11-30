"""
"""

import xml.etree.ElementTree as ET

import requests


OBC_IDE_URL = "http://api.obc.io/publications"


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
        for article_id in root.find("PubmedData").find("ArticleIdList").findall("ArticleId"):
            if article_id.get("IdType") == "pubmed":
                result["pmid"] = article_id.text
            elif article_id.get("IdType") == "doi":
                result["doi"] = article_id.text
            else:
                continue

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

        return result

    def publication_in_obc(self):
        """Return True if an article is found in the obc.ide /publications api response.

        Otherwise, returns False. Checks against pmid first, then doi if pmid is
        not found. Accepts one argument, a PubmedArticle instance.
        """
        entries = requests.get(OBC_IDE_URL).json()

        for entry in entries:
            if self.pmid == entry.get("pmid", None):
                return True
            else:
                if self.doi == entry.get("doi", None):
                    return True
        return False


def sort_xml(path):
    # sort PubmedArticles from specified .xml path into two buckets:
    # those that are found in obc, and those that are not found in obc
    in_obc = []
    not_in_obc = []
    # loop thru the <PubmedArticleSet> node
    tree = ET.parse(path)
    root = tree.getroot()
    for node in root:
        article = PubmedArticle(node)
        if article.publication_in_obc():
            in_obc.append(article)
        else:
            not_in_obc.append(article)
    return in_obc, not_in_obc


if __name__ == '__main__':
    in_obc, not_in_obc = sort_xml('real_output.xml')
    print("In OBC: " + str(len(in_obc)))
    print("Not In OBC: " + str(len(not_in_obc)))
