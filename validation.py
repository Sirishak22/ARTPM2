############################################################
# Given a docs file containing abstracts for results from the 30 topics in topics2017.xml and qrels files and docs,
# does cross-validation of the given LeToR configuration.
#
# Lowell Milliken
############################################################
import random
import learning_to_rank as l2r
import pickle
import os

from learning_to_rank import load_indriscores

cv_dir = 'cv_files'
# m - meta
# sd - splitdrugs
# f - filter
# t - target
# jd - journal disease
# tl - text length
# is - indri scores
features_template = 's_{}_known_features_{}'
unknown_template = 's_{}_{}_unknown_features_{}'
model_name = cv_dir + os.sep + '{}_{}_model_{}'
score_filename = cv_dir + os.sep + '{}_{}_model_{}scores'


def gen_cv_sets():
    cv_sets = [1]*3 + [2]*3 + [3]*3 + [4]*3 + [5]*3 + [6]*3 + [7]*3 + [8]*3 + [9]*3 + [10]*3
    random.seed()
    random.shuffle(cv_sets)
    return cv_sets


# ListNet rparams = {'-lr': 0.1, '-epoch': 3000}
def do_cv(unknown_docs_filename='topics2017_m_as_ex_tr_nd_ft_nsh_prf-2-20-0.5-0.5_large_gfix_alldocs.txt',
          meta=False, splitdrugs=False, metric='P@10', program='RankLib', filtered=False, targetproxy=False, dist=5,
          journaldisease=False, textlen=True, indriscore=False, otherscore=False, ranker='ListNet', rparams=None,
          kscorefile='topics2017_m_as_ex_tr_nd_ft_nsh_prf-2-10-0.5-0.8_basescores_large_gfix_run.txt',
          scorefile='topics2017_m_as_ex_tr_nd_ft_nsh_prf-2-10-0.5-0.8_ob-topics2017_m_as_ex_tr_nd_ft_nsh_prf-2-20-0.5-0.5_large_gfix_large_gfix_run.txt',
          fixparts=True, normscores=False, phraseterms=False, intval=True, termfile=None, termkeyfile=None, nodrugs=False):
    """Creates a ranked list output file in TREC format doing training and cross validation for LeToR.

    :param unknown_docs_filename: name of file containing abstracts from the current Retrieval stage run
    :param meta: Boolean. Use metamap CUIs or not. Requires unknown_docs_filename + '.meta' file containing CUIs for each abstract.
    :param splitdrugs: split drugs into multiple features?
    :param metric: metric to train on. See RankLib help for options.
    :param program: Program to do LeToR with. Default RankLib.
    :param filtered: Filter CUIs to use with meta option. Requires either fterms.pickle or terms_filtered.pickle (for phraseterms) file.
    :param targetproxy: Use proximity to the work 'target' as a feature.
    :param dist: distance threshold for 'target' proximity
    :param journaldisease: Use disease presence in journal name as a feature.
    :param textlen: Use abtract length as a feature.
    :param indriscore: Use the indri score as a feature. Requires Indri scores for the qrel documents called unknown_docs_filename[:-11] + basescores_run.txt and a Indri results file called unknown_docs_filename[:-11] + run.txt
    :param otherscore: Use tf-idf and bm25 scores as a feature. Requires 'qrel_tfidfbase_run.txt' and 'qrel_bm25base_run.txt' as well as unknown_docs_filename[:-11] + 'tfidfbase_run.txt' and unknown_docs_filename[:-11] + 'bm25base_run.txt'
    :param ranker: LeToR ranker to use. See RankLib help for options.
    :param rparams: LeToR parameters in a dictionary. See RankLib help for options. Parameter name including leading '-' is key and parameter value is value.
    :param kscorefile: Alternate score file for use as a feature. This should be scores for the known qrels for training.
    :param scorefile: Alternate score file for use as a feature. This should be scores for the unknown documents for testing.
    :param fixparts: Boolean. Fixed cross-valiation partitions if True.
    :param normscores: Boolean. If True, Indri scores are normalized by (score - minscore)/(maxscore - minscore). Using the '-norm' in rparams with a norm type is preferred. See RankLib help.
    :param phraseterms: Boolean. Use only metamapped CUI terms from original terms that are not unigrams.
    :param intval: Boolean. Use RankLib internal validation. True preferred.
    :param termfile: Explicit set of CUI terms to use. A list in a pickle file.
    :param termkeyfile: Keys for mapping terms in the term file to features. Dict in a pickle file. Key = term. Value = term number (which maps to a feature number).
    :param nodrugs: Boolean. If True, do not use any drug information as a feature.
    """
    unknown_base = unknown_docs_filename[:-11]
    parastr = 'n'

    if meta:
        parastr += '_m'
    if splitdrugs:
        parastr += '_sd'
    if nodrugs:
        parastr += '_nd'
    if filtered:
        parastr += '_f'
    if targetproxy:
        parastr += '_t{}'.format(dist)
    if journaldisease:
        parastr += '_jd'
    if textlen:
        parastr += '_tl'
    if indriscore:
        parastr += '_is'
    if otherscore:
        parastr += '_os'
    if normscores:
        parastr += '_ns'

    if scorefile:
        parastr += '_sf'
    if phraseterms:
        parastr += '_pt'
    if not intval:
        parastr += '_nov'

    topics = l2r.load_topics(distance=dist)

    filteredstr = '_filtered'
    if termfile is None:
        if not phraseterms:
            if filtered:
                termfile = 'terms{}.pickle'.format(filteredstr)
                termkeyfile = 'term_keys{}.pickle'.format(filteredstr)
            else:
                termfile = 'terms{}.pickle'.format('')
                termkeyfile = 'term_keys{}.pickle'.format('')
        else:
            termfile = 'fterms.pickle'
            termkeyfile = 'fterms_keys.pickle'
    else:
        parastr += '_' + termfile[:-11]

    if not os.path.exists(termfile):
        if not filtered:
            meta_docs = l2r.load_docs('qrel_docs.txt.meta')
        else:
            meta_docs = l2r.load_docs('qrel_docs.txt.meta.filtered5')

        l2r.save_terms(meta_docs, filtered)

    with open(termfile, 'rb') as infile:
        terms = pickle.load(infile)
    with open(termkeyfile, 'rb') as infile:
        term_keys = pickle.load(infile)

    meta_docs = None
    unknown_meta_docs = None

    if indriscore:
        basescores = load_indriscores(unknown_base + 'basescores_run.txt', normscores)
        unknownscores = load_indriscores(unknown_base + 'run.txt', normscores)
    else:
        basescores = None
        unknownscores = None

    if otherscore:
        basetfidfscores = load_indriscores('qrel_tfidfbase_run.txt', normscores)
        basebm25scores = load_indriscores('qrel_bm25base_run.txt', normscores)

        unknownitftdfscores = load_indriscores(unknown_base + 'tfidfbase_run.txt', normscores)
        unknownbm25scores = load_indriscores(unknown_base + 'bm25base_run.txt', normscores)
    else:
        basetfidfscores = None
        basebm25scores = None

        unknownitftdfscores = None
        unknownbm25scores = None

    if scorefile:
        kprecscores = load_indriscores(kscorefile, normscores)
        precscores = load_indriscores(scorefile, normscores)
    else:
        kprecscores = None
        precscores = None

    train_all = cv_dir + os.sep + features_template.format(parastr, 'all')
    test_all = cv_dir + os.sep + unknown_template.format(unknown_base, parastr, 'all')
    if filtered:
        train_all += filteredstr
        test_all += filteredstr
    # if not os.path.exists(train_all) or indriscore:
    known_docs = l2r.load_docs()
    if meta:
        if not filtered:
            meta_docs = l2r.load_docs('qrel_docs.txt.meta')
        else:
            meta_docs = l2r.load_docs('qrel_docs.txt.meta.filtered5')

    l2r.save_all_features(topics, known_docs, train_all, known=True, metadocs=meta_docs, terms=terms,
                          term_keys=term_keys, splitdrugs=splitdrugs, targetproxy=targetproxy, journaldisease=journaldisease,
                          textlen=textlen, scores=basescores, tfidfscores=basetfidfscores, bm25scores=basebm25scores,
                          precscores=kprecscores, nodrugs=nodrugs)
    # if not os.path.exists(test_all):
    unknown_docs = l2r.load_docs(unknown_docs_filename)
    if meta:
        unknown_meta_docs = l2r.load_docs(unknown_docs_filename + '.meta')
    l2r.save_all_features(topics, unknown_docs, test_all, known=False, metadocs=unknown_meta_docs, terms=terms,
                          term_keys=term_keys, splitdrugs=splitdrugs, targetproxy=targetproxy, journaldisease=journaldisease,
                          textlen=textlen, scores=unknownscores, tfidfscores=unknownitftdfscores, bm25scores=unknownbm25scores,
                          precscores=precscores, nodrugs=nodrugs)

    cv_file = cv_dir + os.sep + 'cv_sets.txt'
    if fixparts and os.path.exists(cv_file):
        cv_sets = []
        with open(cv_file, 'r') as cvsetfile:
            for line in cvsetfile:
                cv_sets.append(int(line.strip()))
    else:
        cv_sets = gen_cv_sets()
        with open(cv_file, 'w') as cvsetfile:
            for i in cv_sets:
                cvsetfile.write('{}\n'.format(i))

    all_qnos = list(range(1, 31))
    qscores ={}
    pmids = {}
    for i in range(1, 11):
        model_file = model_name.format(parastr, ranker, i)
        train_filename = cv_dir + os.sep + features_template.format(parastr, i)
        test_filename = cv_dir + os.sep + unknown_template.format(unknown_base, parastr, i)
        training_set = [str(x) for x in all_qnos if cv_sets[x-1] != i]
        test_set = [str(x) for x in all_qnos if cv_sets[x-1] == i]

        filter_file(train_all, train_filename, training_set)
        filter_file(test_all, test_filename, test_set)

        # if not os.path.exists(model_file) or indriscore:
        l2r.train_model(train_filename, model_file, ranker=l2r.rankers[ranker], metric=metric, program=program, params=rparams, validation=intval)

        l2r.predict(model_file, test_filename, score_filename.format(parastr, ranker, i), metric=metric, program=program, params=rparams)

        if program == 'RankLib':
            qscores.update(l2r.load_rankings(score_filename.format(parastr, ranker, i)))
            pmids.update(l2r.load_pmids_from_features(test_filename))
        elif program == 'Quickrank':
            qpmids = l2r.load_pmids_from_features(test_filename)
            qscores.update(l2r.load_quickrank_scores(qpmids, score_filename.format(parastr, ranker, i)))
            pmids.update(qpmids)

    runfilename = unknown_base + 'tvs_L2R_{}_{}_{}_run.txt'.format(ranker, metric, parastr)
    l2r.save_reranked(qscores, pmids, runfilename)

    return runfilename


# create new file with only docs
def filter_file(infilename, outfilename, filter_):
    count = 0
    print('')
    with open(infilename, 'r') as infile, open(outfilename, 'w') as outfile:
        for line in infile:
            count += 1
            print('\rFiltering on {}'.format(count), end='')
            qno = line.split()[1].split(':')[1]
            if qno in filter_:
                outfile.write(line)
