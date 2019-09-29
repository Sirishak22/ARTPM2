#######################################
# Converts PubMed results file into TREC format
# Lowell Milliken
#######################################
import os


def main(directory):
    file_base = directory + os.sep + 'pubmed_result({}).txt'
    line_base = '{} Q0 {} {} {} pubmedResult\n'

    with open('pubmed_results.txt', 'w') as outfile:
        for x in range(1, 31):
            with open(file_base.format(x), 'r') as infile:
                results = []
                for line in infile:
                    results.append(line.strip())

                n = len(results)
                for i, pmid in enumerate(results):
                    outfile.write(line_base.format(x, pmid, i+1, n-i))