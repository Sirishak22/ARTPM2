#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 13:35:11 2019

@author: lowellmilliken
"""

abstract="Delineation of the intimate 38 38 38 Liposarcoma of the CDK4 backbone conformation CDK4 GERD GERD GERD GERD of pyridine nucleotide coenzymes in aqueous solution.Fourier Analysis,Magnetic Resonance Spectroscopy,Models, Molecular,Molecular Conformation,NAD,analogs & derivatives,NADP,Structure-Activity Relationship,Temperature"
query="Liposarcoma, CDK4 Amplification, 38-year-old male, GERD"
docset = [query,abstract]
"""
from sklearn.feature_extraction.text import TfidfVectorizer
import pandas as pd

# Create the Document Term Matrix
#count_vectorizer = TfidfVectorizer(stop_words='english')
tfidf_vectorizer = TfidfVectorizer()
sparse_matrix = tfidf_vectorizer.fit_transform(docset)
doc_term_matrix = sparse_matrix.todense()
df = pd.DataFrame(doc_term_matrix, 
                  columns=tfidf_vectorizer.get_feature_names(), 
                  index=['query','abstract'])
print(df)
from sklearn.metrics.pairwise import cosine_similarity
print(cosine_similarity(df, df))
"""
from sklearn.feature_extraction.text import CountVectorizer
import pandas as pd

# Create the Document Term Matrix
#count_vectorizer = TfidfVectorizer(stop_words='english')
count_vectorizer = CountVectorizer()
sparse_matrix = count_vectorizer.fit_transform(docset)
doc_term_matrix = sparse_matrix.todense()
df = pd.DataFrame(doc_term_matrix, 
                  columns=count_vectorizer.get_feature_names(), 
                  index=['query','abstract'])
print(df)
from sklearn.metrics.pairwise import cosine_similarity
print(cosine_similarity(df, df)[0][1])
