###################################################
# Retrieves gene information using EUtilities.
# Lowell Milliken
#
###################################################
from Bio import Entrez
from Bio.Entrez.Parser import ValidationError
import xml.etree.ElementTree as ET
import pickle
import atoms_util
import term_util

_email = 'lmillik@mail.sfsu.edu'
_generifsfile = 'generifs.pickle'


def main():
    pass

# Retrieves and save to a pickle file all Gene RIFs for genes in the topics file using EUtilities.
def save_generifs(filename='topics2017.xml'):
    tree = ET.parse(filename)
    root = tree.getroot()
    rifs = {}
    for topic in root:
        number = topic.attrib['number']
        disease = term_util.clean(topic.find('disease').text)
        gene = term_util.clean(topic.find('gene').text)
        rifs[number] = retrieve_generifs(number, disease, gene)

    with open(_generifsfile, 'wb') as outf:
        pickle.dump(rifs, outf)


# Loads Gene RIFs saved to a pickle file
def load_generifs():
    with open(_generifsfile, 'rb') as inf:
        rifs = pickle.load(inf)
    return rifs


# Retieves Gene RIFs using EUtilities and filter for presence of a disease.
def retrieve_generifs(qno='2', disease='Colon cancer', gene='KRAS G13D BRAF V600E'):
    Entrez.email = _email
    names = atoms_util.unique_names(atoms_util.load_atoms()[qno], 'disease')
    #  names = [disease.lower()]
    term = gene + ' AND ' + disease

    handle = Entrez.esearch(db='gene', term=term, sort='relevance')
    record = Entrez.read(handle)
    handle.close()

    handle = Entrez.efetch(db='gene', id=record['IdList'][0], rettype='xml')
    record = Entrez.read(handle)
    handle.close()

    comments = []

    for comment in record[0]['Entrezgene_comments']:
        if comment['Gene-commentary_type'].attributes['value'] == 'generif':
            if 'Gene-commentary_text' in comment.keys():
                for name in names:
                    if term_util.clean(name) in comment['Gene-commentary_text'].lower():
                        comments.append(comment)
                        break

    return comments


# Save PMIDs from Gene RIFs into a file
def generif_pmids(qno_rifs, outfile='rifdocids.txt'):
    with open(outfile, 'w') as outf:
        for qno, rifs in qno_rifs.items():
            outf.write('qno:{}\n'.format(qno))
            for comment in rifs:
                for ref in comment['Gene-commentary_refs']:
                    if 'Pub_pmid' in ref.keys() and 'PubMedId' in ref['Pub_pmid']:
                        outf.write(ref['Pub_pmid']['PubMedId'] + '\n')


# Retrieve aliases of a gene symbol using EUtilities
def retrieve_aliases(gene='KRAS'):
    Entrez.email = _email
    term = gene
    #print(term)
    handle = Entrez.esearch(db='gene', term=term, sort='relevance')
    record = Entrez.read(handle)
    handle.close()

    handle = Entrez.efetch(db='gene', id=record['IdList'][0], rettype='xml')
    record = Entrez.read(handle)
    handle.close()

    aliases = record[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_syn']
    #print(aliases)
    print(record[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_syn'])
    aliases.append(record[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_desc'])
    return aliases


# Save gene aliases into a pickle file
def save_gene_aliases(genes, genesynfile='genesyn.pickle'):
    gene_syns = {}
    for gene in genes:
        if gene!='':
         try:
            syns = retrieve_aliases(gene)
            gene_syns[gene] = syns
         except ValidationError:
            gene_syns[gene] = []
            print(gene + ' xml parse failed\n')
         except KeyError:
            gene_syns[gene] = []
            print(gene + ' key error\n')

    with open(genesynfile, 'wb') as outfile:
        pickle.dump(gene_syns, outfile)


# load gene aliases from a pickle file.
def load_gene_aliases(genesynfile='genesyn.pickle'):
    with open(genesynfile, 'rb') as infile:
        return pickle.load(infile)


if __name__ == '__main__':
    main()
