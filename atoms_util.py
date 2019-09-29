###########################################################
# Retrieves synonyms (atoms) from UMLS.
# Requires a MetaMap Lite installation. The paths in the mm_lite list must be changed to match the locations in the local install of MetaMap Lite.
# Requires a UMLS API key in the file referred to in the apikeyfile variable.
# Lowell Milliken
###########################################################
import Authentication
import simplejson
import requests
import xml.etree.ElementTree as ET
import subprocess
import pickle
import sys
from term_util import *
import json

apikeyfile = 'umls_api_key'
base_uri = 'https://uts-ws.nlm.nih.gov/rest'
all_atoms = '/content/current/CUI/{}/atoms'
language = 'ENG'
mm_lite = ['java', '-cp', '/Users/lowellmilliken/Documents/precision_medicine_contd/public_mm_lite/target/metamaplite-3.6.2rc3-standalone.jar',
           'gov.nih.nlm.nls.ner.MetaMapLite', '--indexdir=/Users/lowellmilliken/Documents/precision_medicine_contd/public_mm_lite/data/ivf/2018AB/Base/strict',
           '--modelsdir=/Users/lowellmilliken/Documents/precision_medicine_contd/public_mm_lite/data/models', '--specialtermsfile=/Users/lowellmilliken/Documents/precision_medicine_contd/public_mm_lite/data/specialterms.txt']


# Get atoms for CUIs from UMLS.
# tgt is the ticket granting ticket for the UMLS service. 
# TODO: replace api call with local database lookup
def get_atoms(cuis, tgt=None):
    auth = Authentication.Authentication(apikeyfile)

    if not tgt:
        tgt = auth.get_ticket_granting_ticket()
    atoms = {}
    params = {'language': language}
    for cui in cuis:
        #print("In CUI loop")
        #print(cui)
        params['ticket'] = auth.get_service_ticket(tgt)
        request = requests.get(base_uri + all_atoms.format(cui), params=params)
        if request.status_code!=404:
        #dataform = str(request).strip("'<>() ").replace('\'', '\"')
        #results = json.loads(dataform)
        #request = requests.get(base_uri + all_atoms.format(cui), params=params).json()
           #print(request.text)
           results = simplejson.loads(request.text)
        #print(results)

           atoms[cui] = []
           for atom in results['result']:
               atoms[cui].append(atom['name'])
        else:
           atoms[cui]=[]
      

    return atoms


# Get atoms associated with a disease from UMLS.
def get_atoms_from_str(disease):
    auth = Authentication.Authentication(apikeyfile)
    tgt = auth.get_ticket_granting_ticket()
    fn = 'tempfile'
    fnout = fn + '.txt'
    fnin = fn + '.mmi'
    with open(fnout, 'w') as fout:
        fout.write(disease)

    subprocess.run(mm_lite + [fnout])

    cuis = []
    with open(fnin, 'r') as fin:
        for line in fin:
            tokens = line.strip().split('|')
            cuis.append(tokens[4])
    #print("In getatoms method")
    #print(cuis)        
    atoms = get_atoms(cuis, tgt)
    return atoms


# Get all atoms for concepts in a topics XML file and save them in a pickle file.
def save_atoms(filename, atomfile='split_atoms.pickle'):
    tree = ET.parse(filename)
    root = tree.getroot()

    auth = Authentication.Authentication(apikeyfile)
    tgt = auth.get_ticket_granting_ticket()

    qno_atoms = {}

    for topic in root:
        number = topic.attrib['number']
        disease = clean(topic.find('disease').text)
        gene = clean(topic.find('gene').text)
        demographic = clean(topic.find('demographic').text)
        if topic.find('other'):
            other = clean(topic.find('other').text)
        else:
            other = ''

        text = {'disease': disease, 'gene': gene, 'demographic': demographic, 'other': other}
        qno_atoms[number] = {}

        fn = 'tempfile'
        fnout = fn + '.txt'
        fnin = fn + '.mmi'

        for key, value in text.items():
            with open(fnout, 'w') as fout:
                fout.write(value)

            subprocess.run(mm_lite + [fnout])

            cuis = []
            with open(fnin, 'r') as fin:
                for line in fin:
                    tokens = line.strip().split('|')
                    cuis.append(tokens[4])

            atoms = get_atoms(cuis, tgt)
            qno_atoms[number][key] = atoms

    with open(atomfile, 'wb') as outfile:
        pickle.dump(qno_atoms, outfile)


# Load atoms from a pickle file.
def load_atoms(atomfile='split_atoms.pickle'):
    with open(atomfile, 'rb') as infile:
        qno_atoms = pickle.load(infile)

    return qno_atoms

# Get unique names from atoms belonging to the field given or the whole set of atoms if no field is given.
def unique_names(atoms, fieldtopic= None):
    names = set()
    if fieldtopic is not None:
        fatoms = atoms[fieldtopic]
    else:
        fatoms = atoms
    for cui, catoms in fatoms.items():
        for name in catoms:
            if name.lower() not in names:
                names.add(name.lower())
    return names


if __name__ == '__main__':
    save_atoms(sys.argv[1], sys.argv[2])
