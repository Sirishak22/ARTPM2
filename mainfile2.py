#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 25 09:15:11 2019

@author: lowellmilliken
"""

import parse_results_for_top_N
import get_docs_by_file
parse_results_for_top_N.all_pmids_in_file('/Users/lowellmilliken/Downloads/indri-5.14/runquery/output_run.txt','pmidfile.txt')
#get_docs_by_file.main('/Users/lowellmilliken/Downloads/indri-5.14/runquery/output_run.txt', 'getdocsoutput.txt', '/Users/lowellmilliken/Documents/precision_medicine_contd/indexes/medline-ja2018-index')
#pmids = get_docs_by_file.load_pmids('pmidfile.txt')
get_docs_by_file.main('/Users/lowellmilliken/Documents/precision_medicine_contd/lmillik-artpm-c576ced69e03/sampleqrelsfortest.txt', '/Users/lowellmilliken/Documents/precision_medicine_contd/lmillik-artpm-c576ced69e03/doctest2.txt', '/Users/lowellmilliken/Documents/precision_medicine_contd/indexes/medline-ja2018-index-final2')

#print(pmids)