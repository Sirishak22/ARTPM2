###################################################################
# Find the top features in a CoordAsc RankLib 2.10 model file and convert them to CUI terms.
# Lowell Milliken
###################################################################

import sys
import pickle


def main(basefilename, outfileterm, outfilekeys, nfiles=10, filtered=False, offset=10):
    terms = set()
    for n in range(1, nfiles+1):
        terms = terms.union(features_to_terms(basefilename + str(n), filtered=filtered, offset=offset))

    term_keys = {}
    n = 0
    for term in terms:
        term_keys[n] = term
        n += 1

    with open(outfileterm, 'wb') as outfile:
        pickle.dump(terms, outfile)
    with open(outfilekeys, 'wb') as outfile:
        pickle.dump(term_keys, outfile)

    return terms, term_keys


def features_to_terms(modelfile, filtered=False, offset=10):
    filteredstr = '_filtered'
    if filtered:
        termkeyfile = 'term_keys{}.pickle'.format(filteredstr)
    else:
        termkeyfile = 'term_keys{}.pickle'.format('')

    with open(termkeyfile, 'rb') as infile:
        terms = list(pickle.load(infile).keys())

    weights = []
    with open(modelfile, 'r') as infile:
        for line in infile:
            if line.startswith('##'):
                continue

            tokens = line.split(' ')

            for token in tokens:
                fw = token.split(':')
                weights.append((int(fw[0]), float(fw[1])))

    weights.sort(key=lambda tup: abs(tup[1]), reverse=True)

    nterms = []
    minweight = weights[-1][1]
    for fw in weights:
        if fw[1] == minweight:
            break

        feature = fw[0] - offset
        nterms.append(terms[feature])

    return set(nterms)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
