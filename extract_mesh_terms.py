"""
Written in Python 3.5

The purpose of this script is to take a PubMed XML file as input, parse out the MeSH Descriptors and Qualifiers, and count how many of each there is. Output is a tab-delimited text file."""

import xml.etree.ElementTree as ET


def xml_prompt():
    """Prompts user for PubMed XML file."""

    try:
        prompt = input("Please enter the name of the PubMed XML file: ")
        fhand = open(prompt)
        return fhand
    except IOError:
        print("Please ensure that you spelled the file name correctly and that you're in the right directory.")


def check_XML(XML, path):
    """Looks for XML tag and pulls out any text therein."""

    if XML.find(path) is not None:
        if XML.find(path).text is not None:
            return XML.find(path).text
        else:
            return ""
    else:
        return ""


def check_UI(XML, path, attribute):
    """Looks for UI attribute in Descriptor tags and pulls it out."""

    if XML.find(path) is not None:
        for desc in XML.findall(path):
            return desc.get(attribute)
    else:
        return ""


def check_qualifier(XML, path):
    """Looks for XML tag and pulls out the text therein."""

    if XML.find(path) is not None:
        if XML.find(path).text is not None:
            for term in XML.findall(path):
                yield term.text
        else:
            return ""
    else:
        return ""


def parse_descriptor(mesh_tree):
    """Parses MeSH Descriptor Names."""

    descriptor_list = list()

    for term in mesh_tree:
        descriptor_list.append(check_XML(term, './DescriptorName'))
    return descriptor_list


def parse_UI(mesh_tree):
    """Parses Descriptor UIs."""

    UI_list = list()

    for term in mesh_tree:
        UI_list.append(check_UI(term, './DescriptorName', 'UI'))
    return UI_list


def parse_qualifier(mesh_tree):
    """Parses MeSH Qualifier Names."""

    qualifier_list = list()

    for term in mesh_tree:
        qualifier_list.append(check_qualifier(term, './QualifierName'))
    return qualifier_list


def extract_desc(xml_input):
    """Extracts MeSH descriptors, descriptor UIs, and qualifiers from Pubmed DocSum XML."""

    tree = ET.parse(xml_input)
    root = tree.getroot()
    articles = root.findall('./PubmedArticle/MedlineCitation')

    desc_lists = list()
    generator_list = list()
    ui_lists = list()

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
    return generator_list, desc_lists, ui_lists


def count_terms(generators, descriptors, UIs):
    """
    Takes lists of MeSH terms, creates a dictionary to count how many there are, and then prints out the
    results."""

    qual_list = list()
    desc_list = list()
    ui_list = list()

    desc_dict = dict()
    qual_dict = dict()

    for lst in generators:
        deduplicate_list = list()

        for generator in lst:
            for terms in generator:
                duplicate_list = list()

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

    unique_qualifiers = set(qual_list)
    unique_descriptors = set(desc_list)

    print('Number of unique qualifiers: ', len(unique_qualifiers))
    print('Number of unique descriptors: ', len(unique_descriptors))

    dUI_list = [ (zip(ui_list, desc_list)) ]
    dUI_pairs = []

    for lst in dUI_list:
        for element in lst:
            dUI_pairs.append("{}    {}".format(element[0], element[1]))


    # Counts descriptors
    for pair in dUI_pairs:
        desc_dict[pair] = desc_dict.get(pair, 0) + 1

    sorted_desc = [ (value, key) for (key, value) in desc_dict.items() ]

    sorted_desc.sort(reverse = True)

    # Counts qualifiers
    for term in qual_list:
        qual_dict[term] = qual_dict.get(term, 0) + 1

    sorted_qual = [(value, key) for (key, value) in qual_dict.items()]

    sorted_qual.sort(reverse=True)

    formatted_desc = list()
    formatted_qual = list()

    for value, key in sorted_desc:
        formatted_desc.append("%s: %d" % (key, value))

    for value, key in sorted_qual:
        formatted_qual.append("%s: %d" % (key, value))

    return formatted_qual, formatted_desc


if __name__ == '__main__':
    raw = xml_prompt()

    qual, desc, UI = extract_desc(raw)

    qual2, desc2 = count_terms(qual, desc, UI)

    desc_heading = "Desc_UI    Descriptor Counts"
    qual_heading = "Qualifiers Counts"

    with open('mesh_descriptors_test.txt', 'a') as desc_f:
        desc_f.write(desc_heading + '\n\n')
        for pair in desc2:
            desc_f.write("{}\n".format(pair))

    with open('mesh_qualifiers_test.txt', 'a') as qual_f:
        qual_f.write(qual_heading + '\n\n')
        for q in qual2:
            qual_f.write("{}\n".format(q))