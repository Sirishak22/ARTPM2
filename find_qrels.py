# find qrels for specific documents
# usage: python find_qrels.py <document id file>

import sys
import parse_results_for_top_N as pr
import csv


def main(filename, qno='1'):
    qrels = find(filename, qno)
    for qp, qrel in qrels.items():
        print(qp, qrel)


def find(filename, qno='1'):
    qrels = load_qrels()

    target_rels = {}
    with open(filename, 'r') as inFile:
        for line in inFile:
            if line.startswith('qno'):
                qno = line.split(':')[1].strip()
            else:
                pmid = line.strip()
                if pmid in qrels[qno]:
                    target_rels[(qno, pmid)] = qrels[qno][pmid]
                else:
                    target_rels[(qno, pmid)] = 'not listed'

    return target_rels


# returns qrels in [qno:[pmid:relevance]]
def load_qrels(qrelfile='qrels-final-abstracts.txt'):
    qrels = {}
    with open(qrelfile, 'r') as rf:
        for line in rf:
            tokens = line.strip().split()
            if tokens[0] in qrels:
                qrels[tokens[0]][tokens[2]] = tokens[3]
            else:
                qrels[tokens[0]] = {tokens[2]: tokens[3]}
    return qrels


def find_missing_qrels(filename):
    qnos = [str(x) for x in range(1, 31)]
    qrels = load_qrels()
    missing = {}
    for qno in qnos:
        pmids = set(pr.parse_pmids(filename, qno, 1000))
        qrel_pmids = set(list(qrels[qno]))
        missing[qno] = qrel_pmids - pmids

    return missing


def save_pmids(all_pmids, filename):
    with open(filename, 'w') as outfile:
        for qno, pmids in all_pmids.items():
            outfile.write('qno:' + qno + '\n')
            for pmid in pmids:
                outfile.write(pmid + '\n')


def sum_qrels(qrels):
    counts = {}
    for qno, pmids in qrels.items():
        zerocount = 0
        onecount = 0
        twocount = 0
        for pmid, qrel in pmids.items():
            if qrel == '0':
                zerocount += 1
            elif qrel == '1':
                onecount += 1
            elif qrel == '2':
                twocount += 1

        counts[qno] = {'qno': qno, '0': zerocount, '1': onecount, '2': twocount}

    with open('qrel_counts.csv', 'w', newline='') as outfile:
        fieldnames = ['qno', '0', '1', '2']
        outfile.write(','.join(fieldnames) + '\n')
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writerows([value for value in counts.values()])


if __name__ == '__main__':
    if len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 2:
        main(sys.argv[1])
