""" Written in Python 3.5

The purpose of this script is to take a PubMed XML file as input, parse out the MeSH Descriptors and Qualifiers, and
count how many of each there is. Output is a tab-delimited text file."""

import xml.etree.ElementTree as ET
import re


# Prompts user for PubMed XML file
def xml_prompt():
    try:
        prompt = input("Please enter the name of the PubMed XML file: ")
        fhand = open(prompt)
        return fhand
    except IOError:
        print("Please ensure that you spelled the file name correctly and that you're in the right directory.")


# Looks for XML tag and pulls out the text therein
def checkXML(XML, path):
    if XML.find(path) is not None:
        if XML.find(path).text is not None:
            return XML.find(path).text
        else:
            return ""
    else:
        return ""


# Looks for UI attribute in Descriptor tags and pulls it out
def checkUI(XML, path, attribute):
    if XML.find(path) is not None:
        for desc in XML.findall(path):
            return desc.get(attribute)
    else:
        return ""


# Looks for XML tag and pulls out the text therein
def checkQualifier(XML, path):
    if XML.find(path) is not None:
        if XML.find(path).text is not None:
            for term in XML.findall(path):
                yield term.text
        else:
            return ""
    else:
        return ""


# Parses out MeSH Descriptor Names
def parse_descriptor(meshTree):
    descriptor_list = []

    for term in meshTree:
        descriptor_list.append(checkXML(term, './DescriptorName'))

    return descriptor_list


# Parses out Descriptor UIs
def parse_UI(meshTree):
    UI_list = []

    for term in meshTree:
        UI_list.append(checkUI(term, './DescriptorName', 'UI'))

    return UI_list


# Parses out MeSH Qualifier Names
def parse_qualifier(meshTree):
    qualifier_list = []

    for term in meshTree:
        qualifier_list.append(checkQualifier(term, './QualifierName'))

    return qualifier_list


# Extracts MeSH terms from input
def extract_desc(xml_input):
    tree = ET.parse(xml_input)
    root = tree.getroot()
    articles = root.findall('./PubmedArticle/MedlineCitation')

    desc_lists = [] # List of each paper's list of descriptors
    generator_list = []
    ui_lists = [] # List of each paper's list of descriptor UIs

    for article in articles:
        desc_lists.append(parse_descriptor(article.findall('./MeshHeadingList/MeshHeading')))
        generator_list.append(parse_qualifier(article.findall('./MeshHeadingList/MeshHeading')))
        ui_lists.append(parse_UI(article.findall('./MeshHeadingList/MeshHeading')))

    counts = 0

    for lst in desc_lists:
        if len(lst) == 0:
            counts = counts + 1
        else : continue

    print("Number of papers without MeSH terms: ", counts)

    return generator_list, desc_lists, ui_lists # Returns a tuple of the descriptors list and qualifiers list


# Takes lists of MeSH terms, creates a dictionary to count how many there are, and then prints out the results
def count_terms(generators, descriptors, UIs):
    qual_list = []
    desc_list = []  # Single list of all the descriptors from each paper
    ui_list = []  # Single list of all of the descriptor UIs from each paper

    desc_dict = dict()
    qual_dict = dict()

    for lst in generators:
        deduplicate_list = []

        for generator in lst:
            for terms in generator:
                duplicate_list = []

                duplicate_list.append(terms)

                for term in duplicate_list:
                    if term not in deduplicate_list:
                        deduplicate_list.append(term)
                        qual_list.append(term)

    for lst in descriptors:
        deduplicate_list2 = []

        for desc in lst:
            duplicate_list2 = []

            duplicate_list2.append(desc)

            for term in duplicate_list2:
                if term not in deduplicate_list2:
                    deduplicate_list2.append(term)
                    desc_list.append(desc)

    for lst in UIs:
        deduplicate_list3 = []

        for ui in lst:
            duplicate_list3 = []

            ui_list.append(ui)

            for id in duplicate_list3:
                if id not in deduplicate_list3:
                    deduplicate_list3.append(id)
                    ui_list.append(id)

    true_q = set(qual_list)
    true_d = set(desc_list)

    print('Number of unique qualifiers: ', len(true_q))
    print('Number of unique descriptors: ', len(true_d))

    dUI_list = [ (zip(ui_list, desc_list)) ]
    dUI_pairs = []

    for lst in dUI_list:
        for element in lst:
            dUI_pairs.append("{}    {}".format(element[0], element[1]))


    ### Counts descriptors ###
    for pair in dUI_pairs:
        desc_dict[pair] = desc_dict.get(pair, 0) + 1

    sorted_desc = [ (value, key) for (key, value) in desc_dict.items() ] # Creates a list of tuples for each k,v pair

    sorted_desc.sort(reverse = True) # Sorts tuples based on the value

    ### Counts qualifiers ###
    for term in qual_list:
        qual_dict[term] = qual_dict.get(term, 0) + 1

    sorted_qual = [(value, key) for (key, value) in qual_dict.items()] # Creates a list of tuples for each k,v pair

    sorted_qual.sort(reverse=True)  # Sorts tuples based on the value

    fancy_desc = [] # Creates empty list for storing sorted descriptor counts
    fancy_qual = []  # Creates empty list for storing sorted qualifier counts

    for value, key in sorted_desc:
        fancy_desc.append("%s: %d" % (key, value))

    for value, key in sorted_qual:
        fancy_qual.append("%s: %d" % (key, value))

    return fancy_qual, fancy_desc


# Executes program
if __name__ == '__main__':
    raw = xml_prompt()

    qual, desc, UI = extract_desc(raw)

    qual2, desc2 = count_terms(qual, desc, UI)

    desc_heading = "Desc_UI    Descriptor Counts"
    qual_heading = "Qualifiers Counts"

    with open('mesh_descriptorstest.txt', 'a') as desc_f:
        desc_f.write(desc_heading + '\n' + '\n')
        for pair in desc2:
            desc_f.write("{}\n".format(pair))

    with open('mesh_qualifierstest.txt', 'a') as qual_f:
        qual_f.write(qual_heading + '\n' + '\n')
        for q in qual2:
            qual_f.write("{}\n".format(q))