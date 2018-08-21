import xml.etree.ElementTree as ET
import json
from extract_mesh_terms import xml_prompt
from search_grants import checkXML


xml_hand = xml_prompt()
json_hand = xml_prompt()


tree = ET.parse(xml_hand)
root = tree.getroot()
articles = root.findall('./PubmedArticle/MedlineCitation')


xml_titles = []
xml_lst = []


for article in articles:
    xml_dict = dict()

    xml_titles.append(article.find('./Article/ArticleTitle').text)

    xml_dict["title"] = article.find('./Article/ArticleTitle').text
    xml_dict["doi"] = checkXML(article, './Article/ELocationID')
    xml_dict["pmid"] = article.find('./PMID').text

    xml_lst.append(xml_dict)


info = json.loads(json_hand.read())

json_titles = []

for paper in info:
    json_title = paper["pubName"]
    json_titles.append(json_title)
    continue


missing_papers1 = []

for title in xml_titles:
    if title not in json_titles:
        missing_papers1.append(title)
    else:
        continue

missing_papers2 = []

for dictionary in xml_lst:
    for paper in missing_papers1:
        if paper == dictionary["title"]:
            missing_papers2.append(dictionary)


print("Difference: " + str(len(missing_papers1)) + " papers")

with open('papers_to_index_manually.txt', 'a') as missing_f:
    for dictionary in missing_papers2:
        missing_f.write("{}\n".format(dictionary))