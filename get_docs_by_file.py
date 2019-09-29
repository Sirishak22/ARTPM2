# retrieve docs from the Indri index given a file with a list of documents.
# documents from each query should be headed by 'qno:<query number>'
# usage: python <pmidfile> <outputfile> <Indri index directory>

import subprocess
import sys
import os


#Code by Lowell
"""
def main(doclist, outfilename='docs.txt', index='medline-ja2018-index'):
    docnos = load_pmids(doclist)
    os.chdir("/Users/lowellmilliken/Downloads/indri-5.14/dumpindex/")
    owd = os.getcwd()
    print(owd)
    with open(outfilename, 'w') as docs:
        for queryno, docsids in docnos.items():
            print('writing qno:{}'.format(queryno))
            docs.write('<QUERY qno=\'{}\'>\n'.format(queryno))
            
            for docno in docsids:
                diprocess = subprocess.Popen(['./dumpindex', index, 'di', 'docno', docno],
                                             stdout=subprocess.PIPE)
                output, err = diprocess.communicate()
                dtprocess = subprocess.Popen(['./dumpindex', index, 'dt', output], stdout=subprocess.PIPE,
                                             universal_newlines=True)
                output, err = dtprocess.communicate()
                docs.write(str(output) + str('\n'))

            docs.write('</QUERY>')

    print('done')
"""

# Code by sirisha to get all the documents for each query at a time 
def main(doclist, outfilename='docs.txt', index='medline-ja2018-index-final2'):
    docnos = load_pmids(doclist)
    os.chdir("/Users/lowellmilliken/Downloads/indri-5.14/dumpindex")
    #owd = os.getcwd()
    #from subprocess import check_output

    #print(owd)
    with open(outfilename, 'w') as docs:
        for queryno, docsids in docnos.items():
            print('writing qno:{}'.format(queryno))
            docs.write('<QUERY qno=\'{}\'>\n'.format(queryno))
            #cmd= ['./dumpindex', index, 'di', 'docno']
            cmd= ['./dumpindex', index, 'doct', 'docno']
            for ids in docsids:
                cmd.append(ids)
            #print(cmd)
            #docidlist = "10048494 10480505"
            #print(docsids)
            #docidlist=" ".join(docsids)
            #i=0
            #for ids in docsids:
            #    if i!=len(docsids)-1:
            #        docsids[i]="'"+ids+"'"+","
            #        i=i+1
            #    else:
            #        docsids[i]="'"+ids+"'"
            #diprocess = subprocess.Popen(['./dumpindex', index, 'di', 'docno', '10048494', '10523382', '10609119'],
            #                                 stdout=subprocess.PIPE,universal_newlines=True)
            #diprocess = subprocess.Popen(['./dumpindex', index, 'di', 'docno', docidlist],
                                         #stdout=subprocess.PIPE,universal_newlines=True)
            #print(diprocess)
            diprocess = subprocess.Popen(cmd,stdout=subprocess.PIPE,universal_newlines=True)
            output, err = diprocess.communicate()
            #out = check_output(['./dumpindex', index, 'di', 'docno', docidlist])
            
            #result = subprocess.run(['./dumpindex', index, 'di', 'docno', '10048494', '10523382', '10609119' ], stdout=subprocess.PIPE)
            #print(result.stdout)
            #print(out)
            #cmddt=['./dumpindex', index, 'dt']
            #for ids1 in output.split():
            #    cmddt.append(ids1)
            #dtprocess = subprocess.Popen(cmddt, stdout=subprocess.PIPE,
            #                                 universal_newlines=True)
            #output, err = dtprocess.communicate()
            #print(output)
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
