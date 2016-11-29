import unittest
import requests
import xml.etree.ElementTree as ET

import obc


class TestObcPubmed(unittest.TestCase):

    PUBMED_XML_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={}&rettype=xml"

    PARSE_TEST_PMID = "27812667"

    # a pmid from a publication that is returned by obc.ide api
    OBC_IDE_TEST_PMID = "26890470"

    def test_parse_doi_pmid(self):

        expected = {"doi": "10.1590/0037-8682-0145-2016", "pmid": "27812667"}

        res = requests.get(
            self.PUBMED_XML_URL.format(self.PARSE_TEST_PMID)
        ).text

        root = ET.tostring(ET.fromstring(res)[0])
        print(root)

        article_obj = obc.PubmedArticle(root)

        self.assertEqual(article_obj.doi, expected["doi"])
        self.assertEqual(article_obj.pmid, expected["pmid"])

    def test_publication_in_obc(self):
        # testing a pmid lookup
        res = requests.get(self.PUBMED_XML_URL.format(self.OBC_IDE_TEST_PMID)).text
        root = ET.tostring(ET.fromstring(res)[0])

        article_obj = obc.PubmedArticle(root)

        self.assertEqual(
            article_obj.publication_in_obc(), True
        )

    def test_sort_xml(self):
        in_obc, not_in_obc = obc.sort_xml('real_output.xml')
        print(len(in_obc))
        print(len(not_in_obc))
