##########################################################################################
# Takes an xml file containing topics in this format:
#
# <topic number="1">
#   <disease>Liposarcoma</disease>
#   <gene>CDK4 Amplification</gene>
#   <demographic>38-year-old male</demographic>
#   <other>GERD</other>
# </topic>
#
# and generates a IndriRunQuery parameter xml file with a variety of formulations.
#
# Lowell Milliken
##########################################################################################
import os
import pickle
import sys
import xml.etree.ElementTree as ET

import learning_to_rank
from get_docs_by_file import load_pmids
from term_util import *


def main(meta=True, filename='topics2018.xml', isexact=True, istreat=True,
         titleonly=False, issyn=True, isband=False, isscoreif=False, isdemo=False, isfiltdemo=True,
         isDrug=False, drugsyn=True, nodiseaseinC=False, geneinscoreif=False, lociandmut=False, expdrugs=False,
         drug_thres=0, target=False, agej=False, diseasej=False, genediseaset=False, genediseasek=False, nogene=False,
         nomut=False, prf=True, prfparams=(2,20,0.5,0.5), pmidfile=None, baseline=False, large=True, otherbase=False,
         indribase=False, noexp=False):
    """Generate a Indri parameter xml file from a topics xml file.

    :param meta: use metamap results for synonyms (default: True)
    :param filename: topics file name (default: topics2017.xml)
    :param isexact: use exact phrases (default: True)
    :param istreat: use treatment words (default: True)
    :param titleonly: search in title field only (default: False)
    :param issyn: contain synonyms in synonym tag (default: True)
    :param isband: use band tag (default: False)
    :param isscoreif: use scoreif tag for disease filtering (default: False)
    :param isdemo: use demographic information (default: False)
    :param isfiltdemo: filter and bin demographic information is isdemo is True (default: True)
    :param isDrug: use drug expansion (default: False)
    :param drugsyn: contain drugs in synonym tag (default: True)
    :param nodiseaseinC: do not use disease in query (default: False)
    :param geneinscoreif: use scoreif tag for gene filtering (default: False)
    :param lociandmut: use loci and mutation information in query (default: False)
    :param expdrugs: use drugs from Mallory et al. 2015 data (default: False)
    :param drug_thres: confidence threshold for expdrugs (default: 0)
    :param target: look for the word target within 5 words of the genes (default: False)
    :param agej: look for age in journal field (default: False)
    :param diseasej: look for disease in journal field (default: False)
    :param genediseaset: look for gene and disease in title field (default: False)
    :param genediseasek: look for gene and disease in keyword field (default: False)
    :param nogene: do not use gene information (default: False)
    :param nomut: do not use mutation terms such as amplification (default: False)
    :param prf: use pseudo-relevance feedback (default True)
    :param prfparams: PRF parameteres (default: (2,20,0.5,0.5))
    :param pmidfile: file name of pmid file to restrict search (default: None)
    :param baseline: create baseline search params using BM25 and tfidf (default: False)
    :param large: get 5000 results instead of 1000 (default: True)
    :param otherbase: if pmidfile contains PMIDs that are not from the qrels file this should be True (default: False)
    :param indribase: create queries with no expansion or structure (default: False)
    :return:
    """
    # from ???? can't find in TREC 2017 PM papers as of 3/8/2018, but it was definitely in one earlier...
    # The papers were revised between now and then and it may have been removed for some reason from a paper
    # treatment_words = 'surgery resistance therapy recurrence treatment targets prognosis malignancy prognostic study survival therapeutical patient outcome'

    #index_name = 'indexes/medline-ja2018-index'
    index_name='/Users/lowellmilliken/Documents/precision_medicine_contd/indexes/medline-ja2018-index-final2'
    pre_base = '<parameters><index>{}</index><runID>testRun</runID><trecFormat>true</trecFormat>\n'.format(index_name)
    post_base = '</parameters>'

    if prf:
        pre_base = '<parameters><index>{}</index><runID>testRun</runID><trecFormat>true</trecFormat>'.format(index_name)
        pre_base += '<fbDocs>{}</fbDocs><fbTerms>{}</fbTerms><fbMu>{}</fbMu><fbOrigWeight>{}</fbOrigWeight>\n'.format(
            prfparams[0], prfparams[1], prfparams[2], prfparams[3])

    if large:
        pre_base += '<count>5000</count>\n'

    outfilename = filename[:-4] + '_'

    if not indribase:
        if meta:
            outfilename += 'm_'
        else:
            outfilename += 'nm'

        if issyn:
            outfilename += 'as_'

        if isscoreif:
            outfilename += 'sf_'

        if isexact:
            outfilename += 'ex_'

        if istreat:
            outfilename += 'tr_'

        if isdemo and isfiltdemo:
            outfilename += 'fd_'
        elif isdemo:
            outfilename += 'd_'
        else:
            outfilename += 'nd_'

        if titleonly:
            outfilename += 't_'
        else:
            outfilename += 'ft_'

        outfilename += 'nsh_'

        if isDrug:
            outfilename += 'd_'
            if drugsyn:
                outfilename += 'ds_'
                if expdrugs:
                    outfilename += 'ed{}_'.format(drug_thres)

        if nodiseaseinC:
            outfilename += 'ndC_'

        if geneinscoreif:
            outfilename += 'gsf_'

        if lociandmut:
            if not nomut:
                if not nogene:
                    outfilename += 'lnm_'
                else:
                    outfilename += 'lnmng_'
            else:
                if not nogene:
                    outfilename += 'l_'
                else:
                    outfilename += 'lng_'

        if target:
            outfilename += 't_'

        if agej:
            outfilename += 'aj_'
        if diseasej:
            outfilename += 'dj_'
        if genediseaset:
            outfilename += 'gdt_'
        if genediseasek:
            outfilename += 'gdk_'

        if noexp:
            outfilename += 'noexp_'

        if prf:
            outfilename += 'prf-{}-{}-{}-{}_'.format(prfparams[0], prfparams[1], prfparams[2], prfparams[3])

        if pmidfile is not None:
            if otherbase:
                outfilename += 'ob-{}_'.format(pmidfile[:-10])
            else:
                outfilename += 'basescores_'
    else:
        outfilename += 'indribase_'

    if large:
        outfilename += 'large_'

    outfilename += 'gfix_'

    if pmidfile is not None:
        if otherbase:
            outfilename += 'ob-{}_'.format(pmidfile[:-10])
        else:
            outfilename += 'basescores_'

    outfilename += 'params.xml'

    if isDrug:
        if expdrugs:
            topicfile = '{}_drugthres-{}.pickle'.format(filename[:-4], drug_thres)
        else:
            topicfile = '{}_pgkb.pickle'.format(filename[:-4])
    else:
        topicfile = '{}.pickle'.format(filename[:-4])

    if os.path.exists(topicfile):
        with open(topicfile, 'rb') as infile:
            topics = pickle.load(infile)
    else:
        topics = learning_to_rank.load_topics(filename, drug_thres=drug_thres, emdrugs=expdrugs)
        with open(topicfile, 'wb') as outfile:
            pickle.dump(topics, outfile)

    qpmids = None
    if pmidfile is not None:
        qpmids = load_pmids(pmidfile)

    if baseline and pmidfile is not None:
        basename = pmidfile[:-10]
        with open(basename + '_tfidfbase_params.xml', 'w') as tfidffile, open(basename + '_bm25base_params.xml', 'w') as bm25file:
            tfidffile.write(pre_base + '<baseline>tfidf</baseline>\n')
            bm25file.write(pre_base + '<baseline>okapi</baseline>\n')
            #for number in range(1, len(topics) + 1):
            topicscrossval=[1, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 29, 30, 31, 32, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 79, 80]
            for number in topicscrossval:
                topic = topics[str(number)]
                query = ' '.join([topic.disease, ' '.join(topic.disease_syn), ' '.join(topic.genes), ' '.join(topic.other)])
                query = form_query(str(number), clean(query), qpmids)
                tfidffile.write(query)
                bm25file.write(query)

            tfidffile.write(post_base)
            bm25file.write(post_base)

        return

    with open(outfilename, 'w') as outFile:
        outFile.write(pre_base)
        tree = ET.parse(filename)
        root = tree.getroot()

        for topic in root:
            number = topic.attrib['number']
            topic = topics[number]

            if not indribase:
                qstring = generate_query(topic, isexact=isexact, meta=meta, istreat=istreat, isdemo=isdemo,
                                         isfiltdemo=isfiltdemo, issyn=issyn, target=target, isDrug=isDrug, drugsyn=drugsyn,
                                         isband=isband, isscoreif=isscoreif, nodiseaseinC=nodiseaseinC, lociandmut=lociandmut,
                                         nomut=nomut, nogene=nogene, titleonly=titleonly, agej=agej, diseasej=diseasej,
                                         genediseaset=genediseaset, genediseasek=genediseasek, geneinscoreif=geneinscoreif,
                                         noexp=noexp)
            else:
                qstring = clean(topic.base)

            outFile.write(form_query(number, qstring, qpmids))

        outFile.write(post_base)

    print('output to ' + outfilename)



def generate_query(topic, isexact=True, meta=True, istreat=True, isdemo=False, isfiltdemo=True, issyn=True,
                   target=False, isDrug=True, drugsyn=True, isband=False, isscoreif=False, nodiseaseinC=False,
                   lociandmut=True, nomut=False, nogene=False, titleonly=False, agej=False, diseasej=False,
                   genediseaset=False, genediseasek=False, geneinscoreif=False, noexp=False):
    """Generate a query string given a topic.

    :param topic: learning_to_rank.Topic object containing topic information
    :param meta: use metamap results for synonyms (default: True)
    :param isexact: use exact phrases (default: True)
    :param istreat: use treatment words (default: True)
    :param titleonly: search in title field only (default: False)
    :param issyn: contain synonyms in synonym tag (default: True)
    :param isband: use band tag (default: False)
    :param isscoreif: use scoreif tag for disease filtering (default: False)
    :param isdemo: use demographic information (default: False)
    :param isfiltdemo: filter and bin demographic information is isdemo is True (default: True)
    :param isDrug: use drug expansion (default: True)
    :param drugsyn: contain drugs in synonym tag (default: True)
    :param nodiseaseinC: do not use disease in query (default: False)
    :param geneinscoreif: use scoreif tag for gene filtering (default: False)
    :param lociandmut: use loci and mutation information in query (default: True)
    :param target: look for the word target within 5 words of the genes (default: False)
    :param agej: look for age in journal field (default: False)
    :param diseasej: look for disease in journal field (default: False)
    :param genediseaset: look for gene and disease in title field (default: False)
    :param genediseasek: look for gene and disease in keyword field (default: False)
    :param nogene: do not use gene information (default: False)
    :param nomut: do not use mutation terms such as amplification (default: False)
    :return: a string representing the query in the Indri query language.
    """
    # from UTD HLTR
    treatment_words = learning_to_rank.Topic.treatment_words
    drugs = ' '.join([exact(clean(drug)) for druglist in topic.drugs for drug in druglist])
    synonyms = {'disease': '', 'gene': '', 'demographic': '', 'other': ''}

    disease = topic.disease
    if isexact:
        disease = exact(disease)

    if meta:
        synonyms['disease'] = ' '.join([exact(clean(dsyn.strip())) for dsyn in topic.disease_syn])
        synonyms['other'] = ' '.join([exact(clean(osyn.strip())) for osyn in topic.other_syn])
    #print(topic.other)
    if len(topic.other)!=0:
       other = topic.other[0]
       if other == 'none':
        other = ''
        synonyms['other'] = ''
    else:    
        other = ''
        synonyms['other'] = ''
    if istreat:
        treat = treatment_words
    else:
        treat = ''

    if not noexp:
        disease += ' ' + synonyms['disease']
        #            gene += ' ' + synonyms['gene']
        #            gene = clean(gene)
        other += ' ' + synonyms['other']

        if not issyn:
            gene = basegene = ' '.join(
                [exact(clean(g)) + ' ' + ' '.join([exact(clean(ga)) for ga in topic.gene_aliases[i]]) for i, g in
                 enumerate(topic.genes)])
        else:
            gene = basegene = ' '.join(
                [syn(exact(clean(g)) + ' ' + ' '.join([exact(clean(ga)) for ga in topic.gene_aliases[i]])) for i, g in
                 enumerate(topic.genes)])

    else:
        if not issyn:
            gene = basegene = ' '.join([exact(clean(g)) for g in topic.genes])
        else:
            gene = basegene = ' '.join([syn(exact(clean(g))) for g in topic.genes])

    demographic = topic.basedemo
    if not isdemo:
        demographic = ''
    elif isfiltdemo:
        demographic = syn(' '.join(topic.age_group)) + ' ' + syn(' '.join(topic.gender))
    else:
        demographic += ' ' + synonyms['demographic']

    if issyn:
        disease = syn(disease)
        #                gene = syn(gene)
        if not isfiltdemo:
            demographic = syn(demographic)
        other = syn(other)
        if treat:
            treat = syn(treat)

    if target:
        gene += ' #uw5(target ' + syn(gene) + ')'
    if isDrug:
        if drugsyn:
            drugs = syn(drugs)
        gene += ' ' + drugs

    if not isscoreif and isband:
        qstring = band(disease, treat) + '\t' + gene + '\t' + demographic + '\t' + other
    else:
        if nodiseaseinC:
            qstring = gene + '\t' + demographic + '\t' + other + '\t' + treat
        else:
            qstring = disease + '\t' + gene + '\t' + demographic + '\t' + other + '\t' + treat

    if lociandmut:
        if not nomut:
            if not nogene:
                if not noexp:
                    qstring += ' ' + ' '.join([uwindow(syn(exact(clean(g)) + ' ' + ' '.join([exact(clean(ga)) for ga in topic.gene_aliases[i]]))
                                                 + ' ' + clean(topic.loci[i]) + ' ' + clean(topic.mutations[i])) for i, g in
                                         enumerate(topic.genes)])
                else:
                    qstring += ' ' + ' '.join([uwindow(exact(clean(g)) + ' ' + clean(topic.loci[i]) + ' ' +
                                                       clean(topic.mutations[i])) for i, g in enumerate(topic.genes)])

            else:
                qstring += ' ' + ' '.join(
                    [uwindow(clean(topic.loci[i]) + ' ' + clean(topic.mutations[i])) for i, g in enumerate(topic.genes) if
                     topic.loci[i] or topic.mutations[i]])
        else:
            if not nogene:
                qstring += ' ' + ' '.join([uwindow(syn(exact(clean(g)) + ' ' + ' '.join([exact(clean(ga)) for ga in topic.gene_aliases[i]]))
                                             + ' ' + clean(topic.loci[i])) for i, g in enumerate(topic.genes)])
            else:
                qstring += ' ' + ' '.join(
                    [clean(topic.loci[i]) for i, g in enumerate(topic.genes)])

    if titleonly:
        qstring = combine_field(qstring, 'title')
    else:
        qstring = combine(qstring)

    if agej:
        qstring = combine(qstring + ' ' + combine_field(' '.join(topic.age_group), 'journal'))

    if diseasej:
        tempdisease = clean(topic.disease + ' ' + ' '.join(topic.disease_syn))
        qstring = combine(qstring + ' ' + combine_field(tempdisease, 'journal'))

    if genediseaset:
        tempdisease = exact(clean(topic.disease + ' ' + ' '.join(topic.disease_syn)))
        qstring = combine(qstring + ' ' + combine_field(tempdisease + ' ' + basegene, 'title'))

    if genediseasek:
        tempdisease = exact(clean(topic.disease + ' ' + ' '.join(topic.disease_syn)))
        qstring = combine(qstring + ' ' + combine_field(tempdisease + ' ' + basegene, 'keywords'))

    if isscoreif:
        if geneinscoreif:
            qstring = scoreif(uwindow(disease + '\t' + basegene), qstring)
        else:
            qstring = scoreif(disease, qstring)

    return qstring


if __name__ == '__main__':
    if len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        main(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        main()
