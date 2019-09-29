#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 15:19:37 2019

@author: lowellmilliken
"""
import pandas
import pickle
header = ['drug','gene','prob']
df3 = pandas.read_table("/Users/lowellmilliken/Documents/precision_medicine_contd/lmillik-artpm-c576ced69e03/stanford_data/emily/genedrug_relationship_100417_sfsu.tsv", sep='\t',header=0, names=header,usecols=[8,9,12])

df31 = df3.loc[(df3['prob'] >= 0.8)]

df31 = df31[['drug','gene']]

df31 = df31.groupby('gene')['drug'].apply(list)
#df31.to_csv("genedrugdict.csv")
#   df31 = pandas.read_csv("genedrugdict.csv", index_col=0)
genedrugdict = df31.T.to_dict()
pickle.dump(genedrugdict, open("genedrugdict.p", "wb"))  # save it into a file named save.p
#newdict = pickle.load(open("genedrugdict.p", "rb"))
#print(newdict)
header = ['disease','gene','prob']
df1 = pandas.read_table("/Users/lowellmilliken/Documents/precision_medicine_contd/lmillik-artpm-c576ced69e03/stanford_data/emily/genedisease_relationship_100417_sfsu.tsv", sep='\t',header=0, names=header,usecols=[8,9,12])
df11 = df1.loc[(df1['prob'] >= 0.8)]
df11 = df11[['disease','gene']]
df2 = df11.groupby('disease')['gene'].apply(list)
genediseasedict = df2.T.to_dict()
pickle.dump(genediseasedict, open("genediseasedict.p", "wb"))
header = ['entity1name','entity1type','entity2name','entity2type','association']
df5 = pandas.read_table("/Users/lowellmilliken/Documents/precision_medicine_contd/lmillik-artpm-c576ced69e03/stanford_data/pharmgkb/relationships/relationships.tsv", sep='\t',header=0, names=header,usecols=[1,2,4,5,7], dtype={'entity1type':str,'entity2type':str},)
df7 = df5.loc[(df5['entity1type'] == 'Gene') & (df5['entity2type'] == 'Gene')]
df8 = df7.groupby('entity1name')['entity2name'].apply(list)
genegenedict = df8.T.to_dict()
pickle.dump(genegenedict, open("genegenedict.p", "wb"))
df9 = df5.loc[(df5['entity1type'] == 'Gene') & (df5['entity2type'] == 'Disease') & (df5['association'] == 'associated')]
df10 = df5.loc[(df5['entity1type'] == 'Gene') & (df5['entity2type'] == 'Chemical') & (df5['association'] == 'associated')]
df9 = df9[['entity1name','entity2name']]
df10 = df10[['entity1name','entity2name']]
df9 = df9.groupby('entity2name')['entity1name'].apply(list)
df10 = df10.groupby('entity1name')['entity2name'].apply(list)
genedrugdict1 = df9.T.to_dict()
genediseasedict1 = df10.T.to_dict()
from collections import defaultdict
from itertools import chain

d1 = defaultdict(list)
for k, v in chain(genedrugdict1.items(),genedrugdict.items()):
    d1[k].extend(v)

d2 = defaultdict(list)
for k, v in chain(genediseasedict1.items(),genediseasedict.items()):
    d2[k].extend(v)

pickle.dump(d1, open("d1.p", "wb"))
pickle.dump(d2, open("d2.p", "wb"))
pickle.dump(genegenedict, open("genegenedict.p", "wb"))