############################################################
# ARtPM learning to rank methods and Topic loading.
#
# Lowell Milliken
############################################################
from term_util import *
import find_qrels
from EUtilities import load_gene_aliases
from EUtilities import retrieve_aliases
from EUtilities import save_gene_aliases
from gene_drug import DrugGraph
from gene_drug import auto_relationships_file
from gene_drug import load_drug_graphs
import atoms_util
from bs4 import BeautifulSoup
import sys
import xml.etree.ElementTree as ET
import subprocess
import pickle
import itertools
import scholarly
from nltk.stem import PorterStemmer
from nltk.corpus import stopwords
from sklearn.feature_extraction.text import CountVectorizer
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
import time

rankers = {'MART': 0, 'RankNet': 1, 'AdaRank': 3, 'Coordinate Ascent': 4, 'LambdaMART': 6, 'ListNet': 7,
           'Random Forests': 8, 'Linear regression': 9}


# def main(docsfile, model):
#     topics = load_topics()
#     docs = load_docs(docsfile)
#
#     filebase = docsfile[:-11] + model[6:-4] + '_'
#
#     featurefile = filebase + 'features.txt'
#     save_all_features(topics, docs, featurefile, False)
#
#     scorefile = filebase + 'scores.txt'
#
#     predict(model, featurefile, scorefile)
#
#     scores = load_rankings(scorefile)
#     pmids = load_pmids_from_features(featurefile)
#     save_reranked(scores, pmids, filebase + 'L2R_run.txt')


# best so far: ListNet, params = {'-lr': '0.01', '-epoch': '3000'}, metric = 'P@10'
def train_model(train_data, model_file, ranker=4, metric='P@10', program='RankLib', params=None, validation=True):
    """Trains a model using RankLib or QuickRank. Edit jar to the location of RankLib.

    :param train_data: filename of training features file
    :param model_file: output filename for model
    :param ranker: Ranker to use, must be from rankers dict. See RankLib help for details.
    :param metric: Metric to use, default is precision @ 10. See RankLib help for details.
    :param program: RankLib or QuickRank
    :param params: parameters for the ranker in a dict such as: {'-lr': '0.01', '-epoch': '3000'}. See RankLib help for details.
    :param validation: Use internal validation or not. Boolean.
    :return:
    """
    jar = '..\\RankLib-2.10.jar'

    if program == 'RankLib':
        args = ['java', '-Xmx8000m', '-jar', jar, '-train', train_data, '-ranker',
                str(ranker), '-metric2t', metric, '-save', model_file, '-gmax', '2']
        if params:
            for key, value in params.items():
                args.append(key)
                args.append(str(value))
        if validation:
            args.append('-tvs')
            args.append('0.1')

        print(args)

        subprocess.run(args)
    elif program == 'Quickrank':
        met, cutoff = metric.split('@')
        args = ['quickrank/quickrank/bin/quicklearn', '--algo', ranker, '--train', train_data, '--train-metric', met,
                '--train-cutoff', cutoff, '--model-out', model_file]
        subprocess.run(args)


def predict(model, features_file, score_file, program='RankLib', metric='P@10', params=None):
    """Generates scores for input documents based on a previously created model. Edit jar to the location of RankLib.

    :param model: input model file
    :param features_file: features file for documents to rank
    :param score_file: output score file
    :param program: RankLib or QuickRank
    :param metric: metric to use. Default: precision @ 10
    :param params: other parameters for RankLib. Just {'-norm': <>} or not. See RankLib help for details.
    :return:
    """
    jar = '..\\RankLib-2.10.jar'

    if program == 'RankLib':
        # jar = '..\\RankLib-2.1-patched.jar'
        args = ['java', '-Xmx8000m', '-jar', jar, '-load', model, '-rank', features_file,
                '-score', score_file, '-gmax', '2']
        if params:
            if '-norm' in params:
                args.append('-norm')
                args.append(params['-norm'])
        subprocess.run(args)
    elif program == 'Quickrank':
        met, cutoff = metric.split('@')
        args = ['quickrank/quickrank/bin/quicklearn',  '--model-in', model, '--test', features_file, '--test-metric', met,
                '--test-cutoff', cutoff, '--scores', score_file]
        subprocess.run(args)


def load_topics(topicfile='topics2017.xml', distance=0, drug_thres=0.9, emdrugs=False, all_drugs=False):
    """Loads topics from a xml file.

    :param topicfile: path or filename of topics file.
    :param distance: target proximity distance.
    :param drug_thres: confidence threshold for drug data.
    :param emdrugs: Use Mallory et al. 2015 drug data if True.
    :return: Dictionary of topics.
    """
    tree = ET.parse(topicfile)
    root = tree.getroot()

    topics = {}

    if topicfile != 'topics2017.xml':
        atomfile = 'split_atoms_2018.pickle'
    else:
        atomfile = 'split_atoms.pickle'

    if not emdrugs:
        drug_source = 'original'
    elif not all_drugs:
        drug_source = 'auto'
    else:
        drug_source = 'both'
    for topic in root:
        number = topic.attrib['number']
        disease = clean(topic.find('disease').text)
        diseasewithoutclean=topic.find('disease').text
        gene = topic.find('gene').text
        demographic = clean(topic.find('demographic').text)
        demowithoutclean=topic.find('demographic').text
        if topic.find('other') is not None:
            other = clean(topic.find('other').text)
            otherwithoutclean=topic.find('other').text
        else:
            other = ''
            otherwithoutclean=''
        topics[number] = Topic(number, disease, gene, demographic, other,diseasewithoutclean,demowithoutclean,otherwithoutclean, max_distance=distance,
                               drug_source=drug_source, drug_thres=drug_thres, atomfile=atomfile)

    return topics


def load_docs(docfile='qrel_docs.txt'):
    """Loads documents from a file.

    :param docfile: filename of documents file.
    :return: documents in dictionary of dictionary. qno -> pmid -> document
    """
    with open(docfile, 'r', encoding='UTF-8') as infile:
        xml = infile.read()

    xml = '<ROOT>' + xml + '</ROOT>'

    soup = BeautifulSoup(xml, 'html.parser')

    docs = {}

    for query in soup('query'):
        qno = query['qno']
        qdocs = {}
        for sdoc in query('doc'):
            doc = {}

            doc['pmid'] = sdoc.find('docno').text
            if sdoc.find('title'):
                doc['title'] = sdoc.find('title').text
            if sdoc.find('text'):
                doc['abstract'] = sdoc.find('text').text
            if sdoc.find('keywords'):
                doc['keywords'] = sdoc.find('keywords').text
            if sdoc.find('journal'):
                doc['journal'] = sdoc.find('journal').text
            if sdoc.find('authors'):
                doc['authors'] = sdoc.find('authors').text
            if sdoc.find('pubtypes'):
                doc['pubtypes'] = sdoc.find('pubtypes').text
            if sdoc.find('pubdate'):
                doc['pubdate'] = sdoc.find('pubdate').text
            if sdoc.find('language'):
                doc['language'] = sdoc.find('language').text
            if sdoc.find('comments'):
                doc['comments'] = sdoc.find('comments').text
            qdocs[doc['pmid']] = doc

        docs[qno] = qdocs

    return docs


def vectorize(doc, field):
    """Creates a vector of terms from a document.

    :param doc: document
    :param field: field restriction for vector
    :return: dictionary of terms to term counts.
    """
    if field == 'all':
        text = doc.get('title', '') + ' ' + doc.get('abstract', '') + ' ' + doc.get('keywords', '')
    else:
        text = doc.get(field, '')

    terms = {}
    for word in text.split():
        if word in terms:
            terms[word] += 1
        else:
            terms[word] = 1

    return terms


def gen_features(topic, docs, outfile, count, total, qrels, known=True, terms=None, term_keys=None, metadocs=None,
                 splitdrugs=True, targetproxy=False, journaldisease=False, textlen=False, scores=None, tfidfscores=None,
                 bm25scores=None, precscores=None, nodrugs=False):
    """Generate features of a topic.

    :param topic: topics to extract features of
    :param docs: Indri query results documents
    :param outfile: output file object
    :param count: current document count
    :param total: total document count
    :param qrels: known relevances
    :param known: True for documents of known relevance
    :param terms: term list for CUI features
    :param term_keys: term feature numbers for CUI features
    :param metadocs: CUI documents for all query results documents
    :param splitdrugs: split drugs into multiple features?
    :param targetproxy: use target term proximity to gene and a feature
    :param journaldisease: use disease in journal field as a feature
    :param textlen: use document length as a feature
    :param scores: Indri result scores to use a feature
    :param tfidfscores: tfidf result scores to use a feature
    :param bm25scores: BM25 result scores to use a feature
    :param precscores: Alternate result scores to use a feature
    :return:
    """

    qno = topic.qno
    #print("tfidfscores")
    #print(tfidfscores)
    #print("bm25scores")
    #print(bm25scores)
    #print(qno)
    for pmid, doc in docs.items():
        if metadocs:
            if pmid not in metadocs[qno]:
                continue
        count += 1
        if metadocs and terms:
            print('\rworking on: {} of {} with {} terms'.format(count, total, len(terms)), end='')
        else:
            print('\rworking on: {} of {}'.format(count, total), end='')

        if scores:
            score = scores[qno][pmid]
        else:
            score = None

        if len(tfidfscores.keys())>0:
            #print(qno)
            #print(pmid)
            #if pmid not in tfidfscores[int(qno)].keys():
             #   continue
            qno1=int(qno)
            #pmid1=int(pmid)
            score1 = tfidfscores[qno1][pmid]
            #print("tfidfscore")
            #print(tfidfscore)
        else:
            score1 = None

        if len(bm25scores.keys())>0:
            #if pmid not in bm25scores[int(qno)].keys():
            #    continue
            score2 = bm25scores[int(qno)][pmid]
            #print("bm25score")
            #print(bm25score)
        else:
            score2 = None

        if precscores:
            if pmid not in precscores[qno]:
                continue
            precscore = precscores[qno][pmid]
        else:
            precscore = None

        if known:
            outfile.write('{} qid:{} '.format(qrels[qno][pmid], qno))
        else:
            outfile.write('{} qid:{} '.format(0, qno))

        if metadocs:
            features, terms, term_keys, feature_count = topic.create_doc_features(doc, known=known,
                                                                                  terms=terms,
                                                                                  term_keys=term_keys,
                                                                                  metadoc=
                                                                                  metadocs[qno][pmid],
                                                                                  splitdrugs=splitdrugs,
                                                                                  targetproxy=targetproxy,
                                                                                  journaldisease=journaldisease,
                                                                                  textlen=textlen,
                                                                                  score=score,
                                                                                  tfidfscore=tfidfscore,
                                                                                  bm25score=bm25score,
                                                                                  precscore=precscore,
                                                                                  nodrugs=nodrugs)
        else:
           start_time = time.time()
           features, feature_count = topic.create_doc_features(doc, known=known,
                                                               splitdrugs=splitdrugs,
                                                               targetproxy=targetproxy,
                                                               journaldisease=journaldisease,
                                                               textlen=textlen,
                                                               score=score,
                                                               tfidfscore=score1,
                                                               bm25score=score2,
                                                               precscore=precscore,
                                                               nodrugs=nodrugs)
           print("--- %s seconds for each document ---" % (time.time() - start_time))
           """
            features, terms, term_keys, feature_count = topic.create_doc_features(doc, known=known,
                                                                                  splitdrugs=splitdrugs,
                                                                                  targetproxy=targetproxy,
                                                                                  journaldisease=journaldisease,
                                                                                  textlen=textlen,
                                                                                  score=score,
                                                                                  tfidfscore=tfidfscore,
                                                                                  bm25score=bm25score,
                                                                                  precscore=precscore,
                                                                                  nodrugs=nodrugs)
           """
        # feat_string = ''
        # for i in range(1, feature_count + 1):
        #     feat_string += '{}:{} '.format(i, features[i])
        feat_string = ' '.join(['{}:{}'.format(i, features[i]) for i in range(1, feature_count + 1)])
        #print(feat_string)
        # for i in range(feature_count + 1, feature_count + len(terms) + 1):
        #     feat_string += '{}:{} '.format(i, features.get(i, 0))
        if metadocs:
            feat_string += ' ' + ' '.join(['{}:{}'.format(i, features.get(i, 0)) for i in
                                           range(feature_count + 1, feature_count + len(terms) + 1)])

        outfile.write(feat_string)
        outfile.write(' # {}:{}\n'.format(qno, pmid))

    return count


def save_all_features(topics, qdocs, featurefile='known_features.txt', known=True, terms=None, term_keys=None,
                      metadocs=None, splitdrugs=True, targetproxy=False, journaldisease=False, textlen=False, scores=None,
                      tfidfscores=None, bm25scores=None, precscores=None, nodrugs=False):
    """Save features of all topics to a file.

    :param topics: topics to extract and save features of
    :param qdocs: Indri query results documents
    :param featurefile: output filename
    :param known: True for documents of known relevance
    :param terms: term list for CUI features
    :param term_keys: term feature numbers for CUI features
    :param metadocs: CUI documents for all query results documents
    :param splitdrugs: split drugs into multiple features?
    :param targetproxy: use target term proximity to gene and a feature
    :param journaldisease: use disease in journal field as a feature
    :param textlen: use document length as a feature
    :param scores: Indri result scores to use a feature
    :param tfidfscores: tfidf result scores to use a feature
    :param bm25scores: BM25 result scores to use a feature
    :param precscores: Alternate result scores to use a feature
    :return:
    """

    qrels = find_qrels.load_qrels(qrelfile='crossvalidationsetqrels.txt')
    if metadocs and not terms:
        terms = set()
        term_keys = {}
    total = 0
    count = 0
    for qno, docs in qdocs.items():
        total += len(docs)
    print('\n')
    with open(featurefile, 'w') as outfile:
        for qno, docs in qdocs.items():
            count = gen_features(topics[qno], docs, outfile, count, total, qrels, known=known, terms=terms, term_keys=term_keys,
                         metadocs=metadocs, splitdrugs=splitdrugs, targetproxy=targetproxy,
                         journaldisease=journaldisease, textlen=textlen, scores=scores, tfidfscores=tfidfscores,
                         bm25scores=bm25scores, precscores=precscores, nodrugs=nodrugs)
    print('\n')
    return terms, term_keys


def save_terms(metadocs, filtered=False):
    """Extract and save all terms in a documents file.

    :param metadocs: documents
    :param filtered: filter flag (just changes output filename)
    :return:
    """
    terms = set()
    term_keys = {}
    term_count = 0
    for qno, docs in metadocs.items():
        for doc in docs.values():
            doc_vec = vectorize(doc, 'all')

            for term in doc_vec.keys():
                if term not in terms:
                    terms.add(term)
                    term_keys[term] = term_count
                    term_count += 1

    filteredstr = '_filtered'
    if filtered:
        termfile = 'terms{}.pickle'.format(filteredstr)
        termkeyfile = 'term_keys{}.pickle'.format(filteredstr)
    else:
        termfile = 'terms{}.pickle'.format('')
        termkeyfile = 'term_keys{}.pickle'.format('')
    with open(termfile, 'wb') as outfile:
        pickle.dump(terms, outfile)
    with open(termkeyfile, 'wb') as outfile:
        pickle.dump(term_keys, outfile)


def load_pmids_from_features(featurefile):
    """Get PMIDs from a feature file

    :param featurefile: feature filename
    :return: Dictionary of PMIDs by qno.
    """
    pmids = {}
    with open(featurefile, 'r') as infile:
        for line in infile:
            pos = line.find('#')
            comment = line[pos+1:].strip()
            qno, pmid = comment.split(':')
            if qno not in pmids:
                pmids[qno] = []
            pmids[qno].append(pmid)

    return pmids


def load_rankings(scorefile='score.txt'):
    """Load scores from RankLib score output file.

    :param scorefile: score file filename.
    :return: Dictionary of scores. qno -> list of scores for that qno.
    """
    scores = {}
    with open(scorefile, 'r') as infile:
        for line in infile:
            qno, index, score = line.split()
            if qno not in scores:
                scores[qno] = []
            scores[qno].append(score)

    return scores


def save_reranked(qscores, qpmids, outfilename, threshold=5000):
    """Rerank results list based on new scores.

    :param qscores: new scores by qno. qno -> list of scores for that qno
    :param qpmids: PMIDs by qno. qno -> list of PMIDs
    :param outfilename: output results list filename
    :param threshold: Number of documents per qno in output.
    :return:
    """
    trec_template = '{} Q0 {} {} {} l2r\n'

    with open(outfilename, 'w') as outfile:
        for qno in range(1, 31):
            s_pmids = list(zip(qpmids[str(qno)], qscores[str(qno)]))
            s_pmids.sort(key=lambda tup: tup[1], reverse=True)

            s_pmids = s_pmids[:threshold]

            for i, t in enumerate(s_pmids):
                outfile.write(trec_template.format(qno, t[0], i+1, t[1]))


def load_quickrank_scores(qpmids, scorefile):
    """Load scores generated by QuickRank

    :param qpmids: PMIDs by qno. qno -> list of PMIDs
    :param scorefile: scores filename
    :return:
    """
    qnos = []
    pmids = []
    with open(scorefile, 'r') as infile:
        lines = infile.readlines()

    for qno, pms in qpmids.items():
        for p in pms:
            qnos.append(qno)
            pmids.append(p)

    scores = {}
    for qno, pmid, line in zip(qnos, pmids, lines):
        if qno not in scores:
            scores[qno] = []
        scores[qno].append(line.strip())

    return scores


def term_proximity(tokens, term1, term2):
    """Distance between two term in a list of tokens.

    :param tokens: List of tokens.
    :param term1: first term.
    :param term2: second term.
    :return: distance between the terms.
    """
    term1 = term1.lower()
    term2 = term2.lower()
    if term1 in tokens and term2 in tokens:
        term1_indexes = [index for index, value in enumerate(tokens) if value == term1]
        term2_indexes = [index for index, value in enumerate(tokens) if value == term2]
        distances = [abs(item[0] - item[1]) for item in itertools.product(term1_indexes, term2_indexes)]
        return min(distances)
    else:
        return 10000000


class Topic:
    treatment_words = 'prevention prophylaxis prognosis outcome survival treatment therapy personalized'
    treatment = ['prevention','prophylaxis','prognosis','outcome','survival','treatment','therapy','personalized','medication','prescription','remedy','surgery','cure','healing','therapeutics']
    prevention = ['precaution','precautionary','prophylactic']
    prognosis = ['diagnosis','prognostication']
    mutation_words = ['amplification', 'inactivating', 'fusion transcript', 'fusion', 'deletion', 'loss']
    target_words = [PorterStemmer().stem('target')]
    druggraphp, druggraphm = load_drug_graphs()

    def __init__(self, qno, disease, gene, demo, other, diseasewithoutclean,demowithoutclean,otherwithoutclean,max_distance=0, drug_source='original', drug_thres=0.9,
                 atomfile=''):
        """

        :param qno: Query name (any string)
        :param disease: disease string
        :param gene: genes string in the form 'GENE mutation, GENE mutation'
        :param demo: demographic string ex: '67-year-old male'
        :param other: other condition string ex: 'Depression, Heart Disease'
        :param max_distance: maximum distance for 'target' proximity to genes (ignore this)
        :param drug_source: 'original' uses PharmGKB data, 'auto' uses Mallory 2015 data, or 'both' uses both
        :param drug_thres: confidence threshold for Mallory 2015 data, default 0.9
        :param atomfile: location of preprocessed atomfile (not currently used)
        """
        self.qno = qno

        #TODO: refactor to process as topic comes in
        # loads atoms from preprocessed file
        # atoms = atoms_util.load_atoms(atomfile)
        # names = atoms_util.unique_names(atoms[qno], 'disease')

        # Uses MetaMap Lite to get CUIs and requests atoms from 'https://uts-ws.nlm.nih.gov/rest'
        atoms = atoms_util.get_atoms_from_str(disease)
        names = atoms_util.unique_names(atoms, None)

        self.base = disease + ' ' + gene + ' ' + demo + ' ' + other

        disease_syn = []
        for name in names:
            if len(clean(name)) > 2:
                disease_syn.append(name.lower())
        self.disease = disease.lower()
        self.disease_syn = disease_syn

        self.diseasewords = self.disease.split()
        for word in self.disease_syn:
            self.diseasewords.extend(word.split())

        loci = []
        genes = gene.split(',')
        lgenes = []
        for gene in genes:
            pos1 = gene.find('(')
            if pos1 != -1:
                pos2 = gene.find(')')
                locus = gene[pos1 + 1: pos2].lower()
                lgenes.append(gene[:pos1].strip())
            else:
                lgenes.append(gene.strip())
                locus = ''

            loci.append(locus)

        mutations = []
        mgenes = []
        for gene in lgenes:
            gene = gene.strip()
            if gene.isupper():
                mgenes.append(gene.lower())
                mutations.append('')
            else:
                sgene = gene.split()
                mut = []
                g = ''
                for token in sgene:
                    if token.isupper():
                        g = token
                    else:
                        mut.append(token)

                mgenes.append(g.lower())
                mutations.append(' '.join(mut))

            # found = False
            # for mut in Topic.mutation_words:
            #     if mut in gene.lower():
            #         mgenes.append(gene.lower().replace(mut, '').strip())
            #         mutations.append(mut.lower())
            #         found = True
            #         break
            # if not found:
            #     mgenes.append(gene)
            #     mutations.append('')

        #TODO: refactor to process as topic comes in
        all_gene_aliases = load_gene_aliases()
        gene_aliases = []
        #save_gene_aliases(mgenes)
        #print(mgenes)
        for gene in mgenes:
            # from file
             aliases = [x for x in all_gene_aliases.get(gene, []) if len(x) > 2]
            # EUtilities web request
            #if gene!='':
               #aliases = [x for x in retrieve_aliases(gene) if len(x) > 2]
             gene_aliases.append(aliases)

        self.genes = mgenes
        self.loci = loci
        self.mutations = mutations
        self.gene_aliases = gene_aliases

        self.basedemo = demo
        self.age_group = find_age_group(clean(demo))
        self.gender = find_sex(clean(demo))
        self.other = other.lower().split()
        #self.rawdisease=diseasewithoutclean
        #self.rawdemo=demowithoutclean
        #self.rawother=otherwithoutclean
        self.rawquery=diseasewithoutclean + ' ' + gene + ' ' + demowithoutclean + ' ' + otherwithoutclean
        other_syn = []
        if other:
            # names = atoms_util.unique_names(atoms[qno], 'other')
            atoms = atoms_util.get_atoms_from_str(other)
            names = atoms_util.unique_names(atoms, None)
            for name in names:
                if len(clean(name)) > 2:
                    other_syn.append(name.lower())
        self.other_syn = other_syn

        drugs = []
        for gene in self.genes:
            if drug_source == 'original':
                drugs.append(Topic.druggraphp.find_drugs(gene.upper(), drug_thres))
            elif drug_source == 'auto':
                drugs.append(Topic.druggraphm.find_drugs(gene.upper(), drug_thres))
            else:
                drugs.append(Topic.druggraphp.find_drugs(gene.upper(), drug_thres) +
                             Topic.druggraphm.find_drugs(gene.upper(), drug_thres))

        self.drugs = drugs

        self.stemmer = PorterStemmer()
        self.max_distance = max_distance

    # creates a feature list
    # Topic and Document
    # 1. Gene in title
    # 2. Gene in abstract
    # 3. Gene in keywords
    # 4. Disease in title
    # 5. Disease in abstract
    # 6. Disease in keywords
    # Treatment
    # 7. My list of treatment words in title
    # 8. My list of treatment words in in abstract + keywords
    # 9. Stanford drugs associated with gene in doc
    # n-grams (created with metamap)
    def create_doc_features(self, doc, known, terms=None, term_keys=None, metadoc=None, splitdrugs=True,
                            targetproxy=False, journaldisease=False, textlen=False, score=None, tfidfscore=None,
                            bm25score=None, precscore=None, nodrugs=False):
        features = {}
        #print("In Features")
        #print(doc)
        pmid=doc.get('pmid','')
        print(pmid)
        title = doc.get('title', '').lower()
        abstract = doc.get('abstract', '').lower()
        keywords = doc.get('keywords', '').lower()
        journal = doc.get('journal', '').lower()
        authors = [a.strip() for a in doc.get('authors', '').lower().split(',')]
        pubdate = doc.get('pubdate','').lower()
        language = doc.get('language','').lower()
        publist = [p.strip() for p in doc.get('pubtypes', '').lower().split('|')]
        numofcomments = [d for d in doc.get('comments','').lower().split('|')]
        title_stemmed = [clean(self.stemmer.stem(term)) for term in title.lower().split()]
        abstract_stemmed = [clean(self.stemmer.stem(term)) for term in abstract.lower().split()]
        keywords_stemmed = [clean(self.stemmer.stem(term)) for term in keywords.lower().split()]
        
        #print("title")
        #print(title)
        """
        print("Genes")
        print(self.genes)
        print("Disease")
        print(self.disease)
        print("Other")
        print(self.other)
        print("Demo")
        print(self.basedemo)
        print("Base")
        print(self.base)
        print("journal")
        print(journal)
        print("authors")
        print(authors)
        print("pubdate")
        print(pubdate)
        print("language")
        print(language)
        print("publist")
        print(publist)
        print("numofcomments")
        print(numofcomments)
        """
        
        abstract_stemmed=" ".join(abstract_stemmed)
        abstract_stemmed_split=abstract_stemmed.split()
        #print(abstract_stemmed_split)
        base_stemmed=[self.stemmer.stem(term) for term in self.base.lower().split()]
        title_stemmed=" ".join(title_stemmed)
        title_stemmed_split=title_stemmed.split()
        keywords_stemmed=" ".join(keywords_stemmed)
        keywords_stemmed_split=keywords_stemmed.split()
        #print("basestemmed")
        #print(base_stemmed)
        start_time = time.time()
        totalwordmatches=0
        for querywords in base_stemmed:
            if querywords in abstract_stemmed_split:
                totalwordmatches = totalwordmatches + abstract_stemmed_split.count(querywords.lower())
        #print("--- %s seconds for feature1 ---" % (time.time() - start_time))
        features[1]=totalwordmatches
        #start_time = time.time()
        titlequerywordmatches = 0
        for querywords in base_stemmed:
            if querywords in title_stemmed_split:
               titlequerywordmatches = titlequerywordmatches+1
        #print("--- %s seconds for feature2 ---" % (time.time() - start_time))
        features[2]=titlequerywordmatches
        #start_time = time.time()
        abstractquerywordmatches = 0
        for querywords in base_stemmed:
            if querywords in abstract_stemmed_split:
               abstractquerywordmatches = abstractquerywordmatches+1
        #print("--- %s seconds for feature3 ---" % (time.time() - start_time))
        features[3]=abstractquerywordmatches
        #start_time = time.time()
        keywordmatches = 0
        for querywords in base_stemmed:
            if querywords in keywords_stemmed_split:
                keywordmatches = keywordmatches+1
        features[4]=keywordmatches
        #print("--- %s seconds for feature4 ---" % (time.time() - start_time))
        #start_time = time.time()
        count=0
        if abstract == title:
            lengthofabstract = 0
        else:
            count = len(abstract.split())
            lengthofabstract = count
        features[5]=lengthofabstract
        #print("--- %s seconds for feature5 ---" % (time.time() - start_time))
        #start_time = time.time()
        lenofdoc = len(title.split())+len(abstract.split())
        features[6]=lenofdoc
        #print("--- %s seconds for feature6 ---" % (time.time() - start_time))
        #start_time = time.time()
        abtreatment=0
        treatmentwords=[self.stemmer.stem(term) for term in Topic.treatment]
        #print("treatmentwords")
        #print(treatmentwords)
        for treatwords in treatmentwords:
            #treatmentwords = self.stemmer.stem(treatmentwords)
            if treatwords in abstract_stemmed_split:
               abtreatment = abtreatment + 1
        titreatment = 0
        for treatwords in treatmentwords:
            #treatmentwords = self.stemmer.stem(treatmentwords)
            if treatwords in title_stemmed_split:
               titreatment = titreatment + 1
        features[7]=abtreatment+titreatment
        #print("--- %s seconds for feature7 ---" % (time.time() - start_time))
        #start_time = time.time()
        abprevention = 0
        preventwords = [self.stemmer.stem(term) for term in Topic.prevention]
        for prevwords in preventwords:
            if prevwords in abstract_stemmed_split:
               abprevention = abprevention + 1
        tiprevention = 0
        for prevwords in preventwords:
            if prevwords in title_stemmed_split:
               tiprevention = tiprevention + 1
        features[8]=abprevention+tiprevention
        #print("--- %s seconds for feature8 ---" % (time.time() - start_time))
        #start_time = time.time()
        prognosiswords = [self.stemmer.stem(term) for term in Topic.prognosis]
        abprognosis = 0
        for prognowords in prognosiswords:
            if prognowords in abstract_stemmed_split:
                abprognosis = abprognosis + 1
        tiprognosis = 0
        for prognowords in prognosiswords:
            if prognowords in title_stemmed_split:
                tiprognosis = tiprognosis + 1
        features[9]=abprognosis+tiprognosis
        #print("--- %s seconds for feature9 ---" % (time.time() - start_time))
        #start_time = time.time()
        tigenewordmatches = 0
        for genewords in self.genes:
            #print(genewords)
            if genewords.lower() in title_stemmed_split:
                tigenewordmatches = tigenewordmatches+1
        features[10]=tigenewordmatches
        #print("--- %s seconds for feature10 ---" % (time.time() - start_time))
        #start_time = time.time()
        abgenewordmatches = 0
        for genewords in self.genes:
            if genewords.lower() in abstract_stemmed_split:
                abgenewordmatches = abgenewordmatches+1
        features[11]=abgenewordmatches
        #print("--- %s seconds for feature11 ---" % (time.time() - start_time))
        #start_time = time.time()
        tidiseasewordmatches = 0
        diseasestemmed=self.stemmer.stem(self.disease)
        diseasestemmed=" ".join(diseasestemmed)
        if diseasestemmed in title_stemmed_split:
                tidiseasewordmatches = tidiseasewordmatches+1
        features[12]=tidiseasewordmatches
        #print("--- %s seconds for feature12 ---" % (time.time() - start_time))
        #start_time = time.time()
        abdiseasewordmatches = 0
        if diseasestemmed in abstract_stemmed_split:
                abdiseasewordmatches = abdiseasewordmatches+1
        features[13]=abdiseasewordmatches
        #print("--- %s seconds for feature13 ---" % (time.time() - start_time))
        #start_time = time.time()
        d1=pickle.load(open("d1.p", "rb"))
        # Feature 22 & 23
        druggenefeature1 = 0
        druggenefeature2 = 0
        druglist=[]
        for genewords in self.genes:
            for key in d1.keys():
                if key.lower()==genewords.lower():
                    geneword=key
                    druglist.extend(d1[geneword])  
        #print("Drug List")
        #print(druglist)
        for value in list(set(druglist)):
              if value.lower() in abstract.lower():
                 druggenefeature1 = druggenefeature1+1
        features[14]=druggenefeature1
        #print("--- %s seconds for feature14 ---" % (time.time() - start_time))
        #start_time = time.time()
        for values in list(set(druglist)):
              if values.lower() in title.lower():
                 druggenefeature2 = druggenefeature2+1
        features[15]=druggenefeature2
        #print("--- %s seconds for feature15 ---" % (time.time() - start_time))
        #start_time = time.time()
        d1=pickle.load(open("d1.p", "rb"))
        d2=pickle.load(open("d2.p", "rb"))
        flag=False
        for key in d2.keys():
            if key.lower()==self.disease.lower():
                #print(key)
                flag=True
                diseaseword=key
        if flag==True:
           geneset = d2[diseaseword]
        else:
           geneset=[]
        #print("Gene Set")
        #print(geneset)
        # Feature 18 and 19
        Diseasegenefeature1 = 0
        Diseasegenefeature2 = 0
        for value in list(set(geneset)):
          if value.lower() in abstract.lower():
              Diseasegenefeature1 = Diseasegenefeature1+1
        features[16]=Diseasegenefeature1
        #print("--- %s seconds for feature16 ---" % (time.time() - start_time))
        #start_time = time.time()
        for values in list(set(geneset)):
          if values.lower() in title.lower():
              Diseasegenefeature2 = Diseasegenefeature2+1
        features[17]=Diseasegenefeature2
        #print("--- %s seconds for feature17 ---" % (time.time() - start_time))
        #start_time = time.time()
        genegenedict = pickle.load(open("genegenedict.p", "rb"))
        # Feature 16 and 17
        genegenefeature1 = 0
        genegenefeature2 = 0
        genelist=[]
        for genewords in self.genes:
            for key in genegenedict.keys():
                if key.lower()==genewords.lower():
                    geneword1=key
                    genelist.extend(genegenedict[geneword1])
        #print("gene list")
        #print(genelist)
        for value in list(set(genelist)):
              if value.lower() in abstract.lower():
                  genegenefeature1 = genegenefeature1+1
        features[18]=genegenefeature1
        #print("--- %s seconds for feature18 ---" % (time.time() - start_time))
        #start_time = time.time()
        for values in list(set(genelist)):
              if values.lower() in title.lower():
                  genegenefeature2 = genegenefeature2+1
       
        features[19]=genegenefeature2
        #print("--- %s seconds for feature19 ---" % (time.time() - start_time))
        #start_time = time.time()
        intersectedlist=list(set(genelist) & set(geneset))
        #print("intersectedlist")
        #print(intersectedlist)
        #intersectedlist = intersection(genelist,geneset)
        # Feature 20 and 21
        intersectgenefeature1=0
        intersectgenefeature2=0
        for value in list(set(intersectedlist)):
          if value.lower() in abstract.lower():
              intersectgenefeature1 = intersectgenefeature1+1
        features[20]=intersectgenefeature1
        #print("--- %s seconds for feature20 ---" % (time.time() - start_time))
        #start_time = time.time()
        for values in list(set(intersectedlist)):
          if values.lower() in title.lower():
              intersectgenefeature2 = intersectgenefeature2+1
        
        features[21]=intersectgenefeature2
        #print("--- %s seconds for feature21 ---" % (time.time() - start_time))
        #start_time = time.time()
        # feature 24 and 25
        finaldruglist = []
        for genes in intersectedlist:
           for key in d1.keys():
              if key.lower()==genes.lower():
                    #print(key)
                 geneskey=key
                 finaldruglist.extend(d1[geneskey])
        #print("Final drug list")
        #print(finaldruglist)
        finaldruggenefeature1 = 0
        finaldruggenefeature2 = 0
        for value in list(set(finaldruglist)):
          if value.lower() in abstract.lower():
              finaldruggenefeature1 = finaldruggenefeature1+1
        features[22]=finaldruggenefeature1
        #print("--- %s seconds for feature22 ---" % (time.time() - start_time))
        #start_time = time.time()
        for values in list(set(finaldruglist)):
          if values.lower() in title.lower():
              finaldruggenefeature2 = finaldruggenefeature2+1
        #print("--- %s seconds for feature23 ---" % (time.time() - start_time))
        #start_time = time.time()
        features[23]=finaldruggenefeature1
        #print(title)
        """
        features[24]=0
        
        search_query = scholarly.search_pubs_query(title)
        #print(search_query)
        # feature 29
        
        val=next(search_query,None)
        #print(val)
        if val is None:
            #print("No search result")
            features[24]=0
        else:
            #print(val.keys())
            val2=getattr(val,'citedby',None)
            if val2 is not None:    
                val1=val.citedby
                features[24]=val1
            else:
                features[24]=0
        """
        #print("--- %s seconds for feature24 ---" % (time.time() - start_time))
        #start_time = time.time()
        # feature 31,32,33
        #abstractdata = abstract.split()
        #abstractwords = [word.lower() for word in abstractdata]
        absdict={}
        pos=1
        for word in abstract_stemmed_split:
            if word not in absdict.keys():
                absdict[word]=[pos]
                pos=pos+1
            else:
                val=[]
                val=absdict[word]
                val.append(pos)
                absdict[word]=val
                pos=pos+1
        #querydata = self.base.split()
        #querywords = [word.lower() for word in querydata]
        poslist=[]
        for queryword in base_stemmed:
            if queryword in absdict.keys():
                poslist.extend(absdict[queryword])       
        poslist.sort()
        #print("poslist")
        #print(poslist)
        p=0
        featurelist=[]
        while p<len(poslist)-1:
            value=poslist[p+1]-poslist[p]
            featurelist.append(value)
            p=p+1 
        #print("Feature List")
        #print(featurelist)
        minproximity=0
        maxproximity=0
        avgproximity=0
        if len(featurelist)!=0:
           minproximity=min(featurelist)
           maxproximity=max(featurelist)
           avgproximity=sum(featurelist)/len(featurelist)
        features[24]=minproximity
        features[25]=maxproximity
        features[26]=avgproximity
        #print("--- %s seconds for feature25,26,27 ---" % (time.time() - start_time))
        #start_time = time.time()
        text=self.base.split()
        #print(text)
        numofsearchterms=len(set(text))
        features[27]=numofsearchterms
        #print("--- %s seconds for feature28 ---" % (time.time() - start_time))
        #start_time = time.time()
        sumofspecialchars=0
        
        for char in self.rawquery:
           #print(char)
           if char.isalnum() or char == " ":
               sumofspecialchars = sumofspecialchars+0
           else:
               sumofspecialchars = sumofspecialchars+1
        features[28]=sumofspecialchars
        #print("--- %s seconds for feature29 ---" % (time.time() - start_time))
        #start_time = time.time()
        stwords = set(stopwords.words('english')) 
        #print(stwords)
        sw=0
        tiswr=0
        titleterms=title_stemmed_split
        if len(titleterms)!=0:
         for tt in titleterms:
            if tt in stwords:
                sw = sw+1
         tiswr= sw/len(titleterms)
        features[29]=tiswr
        #print("--- %s seconds for feature30 ---" % (time.time() - start_time))
        #start_time = time.time()
        sw=0
        abswr=0
        abstractterms=abstract_stemmed_split
        if len(abstractterms)!=0:
         for at in abstractterms:
            if at in stwords:
                sw = sw+1
         abswr= sw/len(abstractterms)
        features[30]=abswr
        #print("--- %s seconds for feature31 ---" % (time.time() - start_time))
        #start_time = time.time()
        #titleterms=title.split()
        #abstractterms=abstract.split()
        #keywordterms=keywords.split(',')
        #qct=titlequerywordmatches/len(titleterms)
        qct=titlequerywordmatches/len(base_stemmed)
        features[31]=qct
        #print("--- %s seconds for feature32 ---" % (time.time() - start_time))
        #start_time = time.time()
        #qca=abstractquerywordmatches/len(abstractterms)
        qca=abstractquerywordmatches/len(base_stemmed)
        features[32]=qca
        #print("--- %s seconds for feature33 ---" % (time.time() - start_time))
        #start_time = time.time()
        #qck=keywordmatches/len(keywordterms)
        qck=keywordmatches/len(base_stemmed)
        
        
        features[33]=qck
        #print("--- %s seconds for feature34 ---" % (time.time() - start_time))
        #start_time = time.time()
        #print(features)
        if pubdate!='' and pubdate.isdigit():
           age=2019-int(pubdate)
        else:
            age=0
        features[34]=age
        #print("--- %s seconds for feature35 ---" % (time.time() - start_time))
        #start_time = time.time()
        publications=['Adaptive Clinical Trial', 'Address', 'Autobiography', 'Bibliography', 'Biography', 'Case Reports', 'Classical Article', 'Clinical Conference', 'Clinical Study', 'Clinical Trial', 'Clinical Trial, Phase I', 'Clinical Trial, Phase II', 'Clinical Trial, Phase III', 'Clinical Trial, Phase IV', 'Clinical Trial Protocol', 'Clinical Trial, Veterinary', 'Collected Works', 'Comparative Study', 'Congress', 'Consensus Development Conference', 'Consensus Development Conference, NIH', 'Controlled Clinical Trial', 'Dataset', 'Dictionary', 'Directory', 'Duplicate Publication', 'Editorial', 'English Abstract', 'Equivalence Trial', 'Evaluation Studies', 'Expression of Concern', 'Festschrift', 'Government Document', 'Guideline', 'Historical Article', 'Interactive Tutorial', 'Interview', 'Introductory Journal Article', 'Journal Article', 'Lecture', 'Legal Case', 'Legislation', 'Letter', 'Meta-Analysis', 'Multicenter Study', 'News', 'Newspaper Article', 'Observational Study', 'Observational Study, Veterinary', 'Overall', 'Patient Education Handout', 'Periodical Index', 'Personal Narrative', 'Portrait', 'Practice Guideline', 'Pragmatic Clinical Trial', 'Publication Components', 'Publication Formats', 'Publication Type Category', 'Randomized Controlled Trial', 'Research Support, American Recovery and Reinvestment Act', 'Research Support, N.I.H., Extramural', 'Research Support, N.I.H., Intramural', "Research Support, Non-U.S. Gov't","Research Support, U.S. Gov't, Non-P.H.S.", "Research Support, U.S. Gov't, P.H.S.", 'Review', 'Scientific Integrity Review', 'Study Characteristics', 'Support of ResearchSystematic Review', 'Technical Report', 'Twin Study', 'Validation Studies', 'Video-Audio Media', 'Webcasts']
        publications=[x.lower() for x in publications]
        #print("length of publications")
        #print(len(publications))
        i=35
        tot=34+len(publications)
        
        while i<=tot:
            features[i]=0
            i=i+1
        #pubfeatures = [0] * len(publications)
        fv=0

        for pubtype in publist:
           #print(pubtype.text)
          if pubtype!='':
           if pubtype =="comparative study":
               fv=1
           if pubtype in publications:
              d = publications.index(pubtype)
              ind=35+d
              features[ind]=1
           #print("pubindex")
           #print(d)
           
           
        features[110]=fv
        #print("--- %s seconds for feature36-112 ---" % (time.time() - start_time))
        #start_time = time.time()
        languages = ["afr","alb","amh","ara","arm","aze","ben","bos","bul","cat","chi","cze","dan","dut","eng","epo","est","fin","fre","geo","ger","gla","gre","heb","hin","hrv","hun","ice","ind","ita","jpn","kin","kor","lat","lav","lit","mac","mal","mao","may","mul","nor","per","pol","por","pus","rum","rus","san","slo","slv","spa","srp","swe","tha","tur","ukr","und","urd","vie","wel"]
        #langfeatures = [0]*61
        j=111
        tot=110+61
        while j<=tot:
            features[j]=0
            j=j+1
        if language!='':    
         p=languages.index(language)
        #print("langindex")
        #print(p)
         ind1=111+p
         features[ind1]=1
        #print("--- %s seconds for feature112-172 ---" % (time.time() - start_time))
        #start_time = time.time()
       
        #print("--- %s seconds for feature173 ---" % (time.time() - start_time))
        #start_time = time.time()
        jaccardcoefftitle=titlequerywordmatches/(len(base_stemmed)+len(title_stemmed_split))
        features[172]=jaccardcoefftitle
        #print("--- %s seconds for feature174 ---" % (time.time() - start_time))
        #start_time = time.time()
        jaccardcoeffabstract=abstractquerywordmatches/(len(base_stemmed)+len(abstract_stemmed_split))
        features[173]=jaccardcoeffabstract
        #print("--- %s seconds for feature175 ---" % (time.time() - start_time))
        #start_time = time.time()
        jaccardcoeffkeywords=keywordmatches/(len(base_stemmed)+len(keywords_stemmed_split))
        #print(len(keywords_stemmed_split))
        #print(len(abstract_stemmed_split))
        #print(len(base_stemmed))
        
        
        features[174]=jaccardcoeffkeywords
        #print("--- %s seconds for feature176 ---" % (time.time() - start_time))
        #start_time = time.time()
        docset = [self.rawquery,abstract]
        count_vectorizer = CountVectorizer()
        sparse_matrix = count_vectorizer.fit_transform(docset)
        doc_term_matrix = sparse_matrix.todense()
        df = pd.DataFrame(doc_term_matrix, 
                          columns=count_vectorizer.get_feature_names(), 
                          index=['query','abstract'])
        features[175]=cosine_similarity(df, df)[0][1]
        #print("--- %s seconds for feature177 ---" % (time.time() - start_time))
        #start_time = time.time()
        docset1 = [self.rawquery,title]
        count_vectorizer = CountVectorizer()
        sparse_matrix = count_vectorizer.fit_transform(docset1)
        doc_term_matrix = sparse_matrix.todense()
        df = pd.DataFrame(doc_term_matrix, 
                          columns=count_vectorizer.get_feature_names(), 
                          index=['query','title'])
        features[176]=cosine_similarity(df, df)[0][1]
        #print("--- %s seconds for feature178 ---" % (time.time() - start_time))
        #start_time = time.time()
        docset2 = [self.rawquery,keywords]
        count_vectorizer = CountVectorizer()
        sparse_matrix = count_vectorizer.fit_transform(docset2)
        doc_term_matrix = sparse_matrix.todense()
        df = pd.DataFrame(doc_term_matrix, 
                          columns=count_vectorizer.get_feature_names(), 
                          index=['query','keywords'])
        features[177]=cosine_similarity(df, df)[0][1]
        #print("--- %s seconds for feature179 ---" % (time.time() - start_time))
        #start_time = time.time()
        feature_count=177

        features[feature_count+1] = tfidfscore
        feature_count = len(features)
        #print("--- %s seconds for feature180 ---" % (time.time() - start_time))
        #start_time = time.time()
        #if bm25score:
        features[feature_count+1] = bm25score
        #print("--- %s seconds for feature181 ---" % (time.time() - start_time))
        #start_time = time.time()
        feature_count = len(features)
        diseaseinkeywords=0
        keywords=keywords.split(',')
        for words in keywords:
            words=words.lower()
        if self.disease.lower() in keywords:
            diseaseinkeywords=diseaseinkeywords+1
        features[feature_count+1]=diseaseinkeywords
        feature_count = len(features)
        ncomments=0
        for comment in numofcomments:
            if comment!='na' or comment!='':
               ncomments=ncomments+1
        features[feature_count+1]=ncomments
        #print("--- %s seconds for feature182 ---" % (time.time() - start_time))
        #start_time = time.time()
        feature_count=len(features)
            
        #print("features")
        #print(features)
        return features,feature_count
        """
        for i, gene in enumerate(self.genes):
        """ 
            
        """
        title_stemmed = [clean(self.stemmer.stem(term)) for term in title.lower().split()]
        abstract_stemmed = [clean(self.stemmer.stem(term)) for term in abstract.lower().split()]
        keywords_stemmed = [clean(self.stemmer.stem(term)) for term in keywords.lower().split()]

        genes_title = 0
        genes_abstract = 0
        genes_keywords = 0
        for i, gene in enumerate(self.genes):
            in_title = False
            in_abs = False
            in_key = False

            near_target_title = False
            near_target_abstract = False

            if gene in title:
                genes_title += 1
                in_title = True
                if targetproxy:
                    near_target_title = term_proximity(title_stemmed, clean(self.stemmer.stem(gene)),
                                                       Topic.target_words[0]) <= self.max_distance
            if gene in abstract:
                genes_abstract += 1
                in_abs = True
                if targetproxy:
                    near_target_abstract = term_proximity(abstract_stemmed, clean(self.stemmer.stem(gene)),
                                                          Topic.target_words[0]) <= self.max_distance
            if gene in keywords:
                genes_keywords += 1
                in_key = True

            for alias in self.gene_aliases[i]:
                if not in_title:
                    if alias.lower() in title:
                        genes_title += 1
                        if targetproxy:
                            near_target_title = term_proximity(title_stemmed, self.stemmer.stem(clean(alias.lower())),
                                                               Topic.target_words[0]) <= self.max_distance
                        break
            for alias in self.gene_aliases[i]:
                if not in_abs:
                    if alias.lower() in abstract:
                        genes_abstract += 1
                        if targetproxy:
                            near_target_abstract = term_proximity(abstract_stemmed, self.stemmer.stem(clean(alias.lower())),
                                                                  Topic.target_words[0]) <= self.max_distance
                        break
            for alias in self.gene_aliases[i]:
                if not in_key:
                    if alias.lower() in keywords:
                        genes_keywords += 1
                        break

        if len(self.genes) > 0:
            features[1] = genes_title / len(self.genes)
            features[2] = genes_abstract / len(self.genes)
            features[3] = genes_keywords / len(self.genes)
            # features[1] = genes_title
            # features[2] = genes_abstract
            # features[3] = genes_keywords
        else:
            features[1] = 0
            features[2] = 0
            features[3] = 0

        if self.disease in title:
            features[4] = 1
        else:
            found = False
            for word in self.disease_syn:
                if word in title:
                    features[4] = 1
                    found = True
                    break
            if not found:
                features[4] = 0
        if self.disease in abstract:
            features[5] = 1
        else:
            found = False
            for word in self.disease_syn:
                if word in abstract:
                    features[5] = 1
                    found = True
                    break
            if not found:
                features[5] = 0
        if self.disease in keywords:
            features[6] = 1
        else:
            found = False
            for word in self.disease_syn:
                if word in keywords:
                    features[6] = 1
                    found = True
                    break
            if not found:
                features[6] = 0

        # found = False
        # for word in Topic.treatment_words.split():
        #     if word in title:
        #         features[7] = 1
        #         found = True
        #         break
        # if not found:
        #     features[7] = 0
        #
        # found = False
        # for word in Topic.treatment_words.split():
        #     if word in abstract or word in keywords:
        #         features[8] = 1
        #         found = True
        #         break
        # if not found:
        #     features[8] = 0
        #
        found = False
        for word in Topic.treatment_words.split():
            word = self.stemmer.stem(word)
            if word in title_stemmed:
                features[7] = 1
                found = True
                break
        if not found:
            features[7] = 0
    
        found = False
        for word in Topic.treatment_words.split():
            word = self.stemmer.stem(word)
            if word in abstract_stemmed or word in keywords_stemmed:
                features[8] = 1
                found = True
                break
        if not found:
            features[8] = 0

        if not nodrugs:
            if splitdrugs:
                found = False
                for drug in self.drugs:
                    for d in drug:
                        if d.lower() in title:
                            features[9] = 1
                            found = True
                            break
                if not found:
                    features[9] = 0

                found = False
                for drug in self.drugs:
                    for d in drug:
                        if d.lower() in abstract:
                            features[10] = 1
                            found = True
                            break
                if not found:
                    features[10] = 0

                found = False
                for drug in self.drugs:
                    for d in drug:
                        if d.lower() in keywords:
                            features[11] = 1
                            found = True
                            break
                if not found:
                    features[11] = 0
            else:
                found = False
                for drug in self.drugs:
                    for d in drug:
                        if d.lower() in title or d.lower() in abstract or d.lower() in keywords:
                            features[9] = 1
                            found = True
                            break
                if not found:
                    features[9] = 0

        feature_count = len(features)
        if targetproxy:
            if len(self.genes) > 0:
                if near_target_title:
                    features[feature_count + 1] = 1
                else:
                    features[feature_count + 1] = 0
                if near_target_abstract:
                    features[feature_count + 2] = 1
                else:
                    features[feature_count + 2] = 0

                feature_count = len(features)

        if journaldisease:
            found = False
            for word in self.diseasewords:
                if word in journal:
                    features[feature_count + 1] = 1
                    found = True
                    break
            if not found:
                features[feature_count + 1] = 0

            feature_count = len(features)

        if textlen:
            if len(abstract.split()) < 10:
                features[feature_count + 1] = 0
            else:
                features[feature_count + 1] = 1
            feature_count = len(features)

        if score:
            features[feature_count + 1] = score
            feature_count = len(features)

        if tfidfscore:
            features[feature_count + 1] = tfidfscore
            feature_count = len(features)

        if bm25score:
            features[feature_count + 1] = bm25score
            feature_count = len(features)

        if precscore:
            features[feature_count + 1] = precscore
            feature_count = len(features)

        if metadoc:
            # feature_count = max([9] + list(term_keys.values()))
            doc_vec = vectorize(metadoc, 'all')
            doc_len = max(sum(doc_vec.values()), 1)

            for term, count in doc_vec.items():
                if term in term_keys:
                    features[term_keys[term] + feature_count + 1] = count/doc_len

        return features, terms, term_keys, feature_count
    
    """
    def gene_string(self):
        all_gene_str = ''
        for gene, aliases in zip(self.genes, self.gene_aliases):
            gene_str = exact(gene) + ' '
            for alias in aliases:
                gene_str += exact(alias) + ' '
            all_gene_str += syn(gene_str.strip()) + ' '

        return all_gene_str.strip()


# if __name__ == '__main__':
#     main(sys.argv[1], sys.argv[2])


def load_indriscores(indrifile, normscores=False):
    qscores = {}
    with open(indrifile, 'r') as infile:
        for line in infile:
            qno, _, pmid, _, score, _ = line.split()
            if qno not in qscores:
                qscores[qno] = {}
            qscores[qno][pmid] = score

    if normscores:
        nscores = {}
        for qno, pmscores in qscores.items():
            nscores[qno] = {}
            maxscore = float(max(pmscores.values()))
            minscore = float(min(pmscores.values()))
            diff = maxscore - minscore
            for pmid, score in pmscores.items():
                nscores[qno][pmid] = str((float(score) - minscore)/diff)

        qscores = nscores

    return qscores