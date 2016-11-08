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

    def _parse_xml(self, xml_str):
        # parse the data we're interested from the xml (passed as a string)
        # return the data as a dictionary
        # function assumes that xml root is a PubMedArticleSet node,
        # containing only one PubMedArticle
        result = {"pmid": None, "doi": None}
        root = ET.fromstring(xml_str)
        for article_id in root[0].find("PubmedData").find("ArticleIdList").findall("ArticleId"):
            if article_id.get("IdType") == "pubmed":
                result["pmid"] = article_id.text
            elif article_id.get("IdType") == "doi":
                result["doi"] = article_id.text
            else:
                continue

        return

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


def publication_in_obc(article):
    """Return True if an article is found in the obc.ide /publications api response.

    Otherwise, returns False. Checks against pmid first, then doi if pmid is
    not found. Accepts one argument, a PubmedArticle instance.
    """
    entries = requests.get(OBC_IDE_URL).json()

    for entry in entries:
        if article.pmid == entry["pmid"]:
            return True
        else:
            if article.doi == entry["doi"]:
                return True
    return False
