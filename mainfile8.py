#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 16:00:35 2019

@author: lowellmilliken
"""

import parse_results_for_top_N
import get_docs_by_file
import xml_to_params
parse_results_for_top_N.all_pmids_in_file('/Users/lowellmilliken/Downloads/indri-5.14/runquery/output_crossvaltrain.txt','crossvaltestpmidfile.txt',qnos=(1, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 29, 30, 31, 32, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 79, 80))
#get_docs_by_file.main('/Users/lowellmilliken/Downloads/indri-5.14/runquery/output_run.txt', 'getdocsoutput.txt', '/Users/lowellmilliken/Documents/precision_medicine_contd/indexes/medline-ja2018-index')
#pmids = get_docs_by_file.load_pmids('pmidfile.txt')
get_docs_by_file.main('/Users/lowellmilliken/Documents/precision_medicine_contd/lmillik-artpm-c576ced69e03/crossvaltestpmidfile.txt', '/Users/lowellmilliken/Documents/precision_medicine_contd/lmillik-artpm-c576ced69e03/crossvalidationtestalldocs.txt', '/Users/lowellmilliken/Documents/precision_medicine_contd/indexes/medline-ja2018-index-final2')
xml_to_params.main(filename='crossvalidationtopics.xml', baseline=True, pmidfile='crossvaltestpmidfile')
