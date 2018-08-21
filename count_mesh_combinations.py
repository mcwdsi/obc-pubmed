import xml.etree.ElementTree as ET


def get_descriptor(mesh_heading, mesh_descriptor, qualifier=None):
    mesh_term = ''

    if mesh_heading.find('./DescriptorName') is not None:
        if mesh_heading.find('./DescriptorName').text is not None:
            desc = mesh_heading.find('./DescriptorName').text

            if desc == mesh_descriptor:
                mesh_term += desc

                if qualifier is not None:
                    qual = get_qualifier(mesh_heading, qualifier)

                    mesh_term += '/'
                    mesh_term += qual
    return mesh_term


def get_qualifier(mesh_heading, mesh_qualifier):
    if mesh_heading.find('./QualifierName') is not None:
        if mesh_heading.find('./QualifierName').text is not None:
            for heading in mesh_heading.findall('./QualifierName'):
                if heading.text == mesh_qualifier:
                    return mesh_qualifier
        else:
            return ""
    else:
        return ""


fhand = input("Location of Pubmed DocSum XML file: ")

tree = ET.parse(fhand)
root = tree.getroot()
articles = root.findall('./PubmedArticle/MedlineCitation')

# Python lists for each MeSH term
comp_sim_lst = list()
software_lst = list()
b_models_lst = list()
t_models_lst = list()
at_least_one = set()

# Python lists for combinations
b_models_and_communicable_disease_lst = list()
b_models_and_infection_lst = list()
t_models_and_communicable_disease_lst = list()
t_models_and_infection_lst = list()
comp_sim_and_communicable_disease_lst = list()
comp_sim_and_infection_lst = list()
software_and_communicable_disease_lst = list()
software_and_infection_lst = list()

for article in articles:
    pmid = article.find('PMID').text
    article_mesh_headings = article.findall('./MeshHeadingList/MeshHeading')

    for heading in article_mesh_headings:

        # Computer Simulation
        if get_descriptor(heading, 'Computer Simulation'):
            comp_sim_lst.append(pmid)
            at_least_one.add(pmid)

        # Software
        if get_descriptor(heading, 'Software'):
            software_lst.append(pmid)
            at_least_one.add(pmid)

        # Biological models
        if get_descriptor(heading, 'Models, Biological'):
            b_models_lst.append(pmid)
            at_least_one.add(pmid)

        # Theoretical models
        if get_descriptor(heading, 'Models, Theoretical'):
            t_models_lst.append(pmid)
            at_least_one.add(pmid)

    # No longer counting how many were assigned 'Communicable Diseases' or 'Infection', in combination
    # with either of the other four descriptors.
    """for heading in article_mesh_headings:
        # Biological models + Communicable diseases / Theoretical models + Communicable diseases /
        # Software + Communicable diseases / Computer simulation + Communicable diseases
        if get_descriptor(heading, 'Communicable Diseases'):
            if pmid in b_models_lst:
                b_models_and_communicable_disease_lst.append(pmid)
            if pmid in t_models_lst:
                t_models_and_communicable_disease_lst.append(pmid)
            if pmid in software_lst:
                software_and_communicable_disease_lst.append(pmid)
            if pmid in comp_sim_lst:
                comp_sim_and_communicable_disease_lst.append(pmid)

        # Biological models + Infection / Theoretical models + Infection /
        # Software + Infection / Computer simulation + Infection
        if get_descriptor(heading, 'Infection'):
            if pmid in b_models_lst:
                b_models_and_infection_lst.append(pmid)
            if pmid in t_models_lst:
                t_models_and_infection_lst.append(pmid)
            if pmid in software_lst:
                software_and_infection_lst.append(pmid)
            if pmid in comp_sim_lst:
                comp_sim_and_infection_lst.append(pmid)"""

print("Software: ", len(software_lst))
print("Computer Simulation: ", len(comp_sim_lst))
print("Models, Theoretical: ", len(t_models_lst))
print("Models, Biological: ", len(b_models_lst))
print("Number of publications with at least one of the descriptors: ", len(at_least_one))
