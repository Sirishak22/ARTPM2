Readme that I sent to Eric for LeToR steps.
All python script ran using 3.6

Software needed:
	RankLib-2.10.jar (https://sourceforge.net/projects/lemur/files/lemur/RankLib-2.10/)
	MetaMap Lite (https://metamap.nlm.nih.gov/MetaMapLite.shtml and need to create an account)

Python modules:
	nltk
	Beautiful Soup (bs4)
	networkx
	BioPython (Bio)

Edit the atom_util.py file:
	mm_lite variable needs to point to location of metamaplite jar file and other file in the directory.
Edit learning_to_rank.py:
	line 86: jar variable should be location of RankLib jar file

Below is the full pipeline from initial query to feature generation to learning to rank.
You can call IndriRunQuery from a python script using subprocess, so you might want to do that to streamline things.
---------------------
Run the Indri query and get a results file with 5000 results per topic: IndriRunQuery command > results file name

In a python script: parse_results_for_top_N.all_pmids_in_file(results file name, pmid file name, number of queries in results file)

Run: python get_docs_by_file.py [pmidfile] [docs file] [Indri index directory]

In a python script:
topic = learning_to_rank.Topic(disease, genes, demo, other) # for each topic
docs = learning_to_rank.load_docs(docs file)[query number] # query number = '1' for first in docs file, '2' for second, etc.
with open('feature file name', 'w') outfile:
	learning_to_rank.save_features(topic, docs, outfile, 0, len(docs), qrels=None, known=False, splitdrugs=False, textlen=True, precscores=None) # precscores will come back later after we get to the query generation component

learning_to_rank.predict('model file name', 'feature file name', 'output score file name')
pmids = learning_to_rank.load_pmids_from_features('feature file name')
scores = learning_to_rank.load_rankings('score file name')
learning_to_rank.save_reranked(scores, pmids, 'final results file name', threshold=1000)
