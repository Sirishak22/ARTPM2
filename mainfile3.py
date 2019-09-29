#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 07:54:29 2019

@author: lowellmilliken
"""

import find_qrels
qrels=find_qrels.load_qrels('/Users/lowellmilliken/Documents/precision_medicine_contd/lmillik-artpm-c576ced69e03/crossvalidationsetqrels.txt')
find_qrels.save_pmids(qrels,'pmidforcrossvalidationsetqrels.txt')