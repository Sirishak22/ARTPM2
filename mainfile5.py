#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 16:05:21 2019

@author: lowellmilliken
"""

import parse_results_for_top_N
tfidfscore=parse_results_for_top_N.load_indri_tfidf_scores(res_filename='/Users/lowellmilliken/Downloads/indri-5.14/runquery/output_tfidf_run.txt')
bm25score=parse_results_for_top_N.load_indri_tfidf_scores(res_filename='/Users/lowellmilliken/Downloads/indri-5.14/runquery/output_bm25_run.txt')
print(tfidfscore)
import time
import learning_to_rank

topics = learning_to_rank.load_topics('topics201802.xml', all_drugs=True)
#print(topics)
#i=1
docs = learning_to_rank.load_docs('doctest2.txt') 
with open('featurefile3.txt', 'w') as outfile:
   for topic in topics.values():
       #print(topic.qno)
       # query number = '1' for first in docs file, '2' for second, etc.
       #print(docs.keys())
       #print(docs['1'])
       #print(docs)
       #with open('featurefile.txt', 'a') as outfile:
          #learning_to_rank.gen_features(topic, docs['1'], outfile, 0, len(docs), qrels=None, known=False, splitdrugs=False, textlen=True, precscores=None) # precscores will come back later after we get to the query generation component
       start_time = time.time()
       learning_to_rank.gen_features(topic, docs[topic.qno], outfile, 0, len(docs), qrels=None, known=False, splitdrugs=False, textlen=True, precscores=None,tfidfscores=tfidfscore,bm25scores=bm25score)
       print("--- %s seconds for each topic ---" % (time.time() - start_time))



