# Modified get_docs_by_file to use multiple indexes
# command line no longer works, use by calling main from a script for interactive python console
# retrieve docs from the Indri index given a file with a list of documents.
# documents from each query should be headed by 'qno:<query number>'
# usage: python <pmidfile> <outputfile> <Indri index directory>

import subprocess
import sys


def main(doclist, outfilename='docs.txt', indexes=['index1', 'index2', 'extraindex'], thresholds=[500, 1000]):
    docnos = load_pmids(doclist)

    with open(outfilename, 'w') as docs:
        for queryno, docsids in docnos.items():
            print('writing qno:{}'.format(queryno))
            docs.write('<QUERY qno=\'{}\'>\n'.format(queryno))
            for docno in docsids:
                if docno.isdigit():
                    for i, threshold in enumerate(thresholds):
                        if int(docno) <= threshold:
                            index = indexes[i]
                            break
                else:
                    index = indexes[len(indexes) - 1]
                diprocess = subprocess.Popen(['dumpindex', index, 'di', 'docno', docno],
                                             stdout=subprocess.PIPE)
                output, err = diprocess.communicate()
                dtprocess = subprocess.Popen(['dumpindex', index, 'dt', output], stdout=subprocess.PIPE,
                                             universal_newlines=True)
                output, err = dtprocess.communicate()
                docs.write(str(output) + str('\n'))

            docs.write('</QUERY>')

    print('done')


def load_pmids(doclist):
    docnos = {}

    with open(doclist, 'r') as docfile:
        for line in docfile:
            if line.startswith('qno'):
                qno = line.split(':')[1].strip()
                docnos[qno] = []
            else:
                docnos[qno].append(line.strip())

    return docnos


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
