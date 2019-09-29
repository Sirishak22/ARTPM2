########################################################################
# Extract doc IDs from TREC format results file.
# Lowell Milliken
#
# Get top N (10) pmids from a query results file for a specific query.
# usage: python parse_results_for_top_N.py <results file> <query name>
########################################################################
import sys


def main(filename, qno):
    pmids = parse_pmids(filename, qno)
    for pmid in pmids:
        print(str(pmid))


# Extract doc IDs for a query from file
def parse_pmids(filename, qno, n=10):
    with open(filename, 'r') as inFile:
        lines = []
        for line in inFile:
            if line.startswith(qno + ' Q'):
                tokens = line.split()
                lines.append((tokens[2], int(tokens[3])))

    lines.sort(key=lambda tup: tup[1])

    if n <= len(lines):
        return [x[0] for x in lines[:n]]
    else:
        return [x[0] for x in lines]


# Extract and save docIDs from a results file into another file
def all_pmids_in_file(res_filename, pmid_filename, qnos=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50)):

    with open(pmid_filename, 'w') as outfile:
        for i in qnos:
            pmids = parse_pmids(res_filename, str(i), 10000)
            outfile.write('qno:{}\n'.format(i))
            for pmid in pmids:
                outfile.write('{}\n'.format(pmid))

def load_indri_tfidf_scores(res_filename):
    tfidfscores={}
    with open(res_filename, 'r') as inFile:
        #lines = []
        for line in inFile:   
            tokens = line.split()
            #tfidfscores[tokens[0]][tokens[2]]=tokens[4]
            
            #lines.append(tokens[0])
            if int(tokens[0]) not in tfidfscores:
                tfidfscores[int(tokens[0])]={}
                #tfidfscores[int(tokens[0]]=tokens[2]
                tfidfscores[int(tokens[0])][tokens[2]]=float(tokens[4])
            else:
                tfidfscores[int(tokens[0])][tokens[2]]=float(tokens[4])
        

    #print(tfidfscores['1'])
    #print(tfidfscores)
    #print(tfidfscores[1])
    #print(tfidfscores[1][10065107])
    return tfidfscores       
    #print(tfidfscores)
    
    
if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
