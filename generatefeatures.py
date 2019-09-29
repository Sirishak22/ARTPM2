#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 07:16:00 2019

@author: lowellmilliken
"""
import learning_to_rank
import time
topics = learning_to_rank.load_topics('topics201801.xml', all_drugs=True)
#print(topics)
#i=1
docs = learning_to_rank.load_docs('doctest1.txt') 
with open('featurefile.txt', 'w') as outfile:
   for topic in topics.values():
       #print(topic.qno)
       # query number = '1' for first in docs file, '2' for second, etc.
       #print(docs.keys())
       #print(docs['1'])
       #print(docs)
       #with open('featurefile.txt', 'a') as outfile:
          #learning_to_rank.gen_features(topic, docs['1'], outfile, 0, len(docs), qrels=None, known=False, splitdrugs=False, textlen=True, precscores=None) # precscores will come back later after we get to the query generation component
    
       start_time = time.time()
       learning_to_rank.gen_features(topic, docs[topic.qno], outfile, 0, len(docs), qrels=None, known=False, splitdrugs=False, textlen=True, precscores=None)
       print("--- %s seconds for each topic ---" % (time.time() - start_time))

