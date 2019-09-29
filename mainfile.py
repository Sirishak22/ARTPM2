#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 13:59:42 2019

@author: lowellmilliken
"""

import learning_to_rank
#from gene_drug import save_all 
#import gene_drug
import xml_to_params
import term_util
topicfile = 'topics201801.xml'

# returns a dict of <topic number> -> <topic object>
#gene_drug.save_drug_graph('drug_data/relationships/relationships.tsv', 'original', 'pharmgkbDG.pickle')
#save_all()
topics = learning_to_rank.load_topics(topicfile, all_drugs=True)
print(topics)

for topic in topics.values():
    
    qstring = xml_to_params.generate_query(topic)
    qstring = term_util.form_query(topic.qno, qstring)
    print(qstring)