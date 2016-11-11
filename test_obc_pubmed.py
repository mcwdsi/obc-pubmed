import unittest
import requests

import obc_pubmed as obc


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

        article_obj = obc.PubmedArticle(res)

        self.assertEqual(article_obj.doi, expected["doi"])
        self.assertEqual(article_obj.pmid, expected["pmid"])

    def test_publication_in_obc(self):
        # testing a pmid lookup
        article_obj = obc.PubmedArticle(
            requests.get(self.PUBMED_XML_URL.format(self.OBC_IDE_TEST_PMID)).text
        )

        self.assertEqual(
            article_obj.publication_in_obc(), True
        )
