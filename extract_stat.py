############################################
# gets metrics from trec statistics file
# Lowell Milliken
############################################
import sys
import subprocess
import os

# Create stats files using trec_eval_9.exe and sample_eval.pl using TREC 2017 Precision Medicine relevance judgments
def create_stats(filename):
    filenamebase = filename[:-8]
    statsfile = filenamebase + '_stats.txt'
    statsfile2 = filenamebase + '_samplestats.txt'

    args = ['trec_eval_9.exe', '-q', '-m', 'all_trec', 'qrels-final-abstracts.txt', filename]

    with open(statsfile, 'w') as statsfilehandle:
        subprocess.run(args, stdout=statsfilehandle)

    args = ['perl', 'sample_eval.pl', '-q', 'sample-qrels-final-abstracts.txt', filename]
    with open(statsfile2, 'w') as statsfilehandle:
        subprocess.run(args, stdout=statsfilehandle)

    return [statsfile, statsfile2]

# Load average from stats file for metrics in "stats"
def get_stats(statsfile, stats):
    avgstats = {}
    with open(statsfile, 'r') as inFile:
        for line in inFile:
            if line.startswith(stats):
                tokens = line.strip().split()
                if tokens[1] == 'all':
                    avgstats[tokens[0]] = tokens[2]

    return avgstats

# Load query level and average stats from stats file for metrics in "stats"
def get_all_stats(statsfile, stats):
    allstats = {}
    for stat in stats:
        allstats[stat] = []

    qnos = []

    with open(statsfile, 'r') as inFile:
        for line in inFile:
            tokens = line.strip().split()
            if tokens[0] in stats:
                allstats[tokens[0]].append((tokens[1], tokens[2]))
                if tokens[1] not in qnos:
                    qnos.append(tokens[1])

    for stat, values in allstats.items():
        values.sort(key=key)

    qnos.sort(key=strkey)
    return allstats, qnos


# put query level stats for all runs in a directory into a csv file.
def all_paired(directory, qfilename='querystats3.csv'):
    files = os.listdir(directory)
    for file in files:
        if file.endswith('_run.txt'):
            paired_stats(directory + os.sep + file, qfilename)


# Get query level stats from a run in "filename" and save to a csv file.
def paired_stats(filename, qfilename='querystats5.csv'):
    filenamebase = filename[:-8]
    stats = ('P_10', 'P_30', 'P_100', 'recall_1000', 'num_rel_ret', 'num_rel', 'ndcg_cut_1000', 'map_cut_1000', 'Rprec_mult_0.20',
             'Rprec_mult_0.40', 'Rprec_mult_0.60', 'Rprec_mult_0.80', 'Rprec_mult_1.00', 'Rprec_mult_1.20',
             'Rprec_mult_1.40', 'Rprec_mult_1.60', 'Rprec_mult_1.80', 'Rprec_mult_2.00', 'infNDCG',
             'iprec_at_recall_0.00', 'iprec_at_recall_0.10', 'iprec_at_recall_0.20', 'iprec_at_recall_0.30',
             'iprec_at_recall_0.40', 'iprec_at_recall_0.50', 'iprec_at_recall_0.60', 'iprec_at_recall_0.70',
             'iprec_at_recall_0.80', 'iprec_at_recall_0.90', 'iprec_at_recall_1.00')

    if not os.path.exists(qfilename):
        with open(qfilename, 'w') as qstatsfile:
            qstatsfile.write('formulation,')
            for qno in range(1, 31):
                qstatsfile.write(','.join(['{}:{}'.format(qno, stat) for stat in stats]) + ',Recall_5000,')
            qstatsfile.write(','.join(['all:{}'.format(stat) for stat in stats]) + ',Recall_5000,')
            qstatsfile.write('\n')

    statsfile, statsfile2 = create_stats(filename)

    allstats, qnos = get_all_stats(statsfile, stats)

    allstats['infNDCG'] = get_all_stats(statsfile2, stats)[0]['infNDCG']

    re5000s = []
    with open(qfilename, 'a+') as qstatsfile:
        qstatsfile.write(os.path.basename(filenamebase))
        for qno in qnos:
            for stat in stats:
                if qno != 'all':
                    qstatsfile.write(',{}'.format(allstats[stat][int(qno) - 1][1]))
                else:
                    qstatsfile.write(',{}'.format(allstats[stat][len(qnos) - 1][1]))
            if qno != 'all':
                re5000 = int(allstats['num_rel_ret'][int(qno) - 1][1]) / int(allstats['num_rel'][int(qno) - 1][1])
                re5000s.append(re5000)
            else:
                re5000 = sum(re5000s) / len(re5000s)
            qstatsfile.write(',{}'.format(re5000))

        qstatsfile.write('\n')

    return

    p10s = []
    num_rel_rets = []
    num_rels = []
    ndcgs = []
    maps = []
    recalls = []
    rprecs = []
    infndcgs = []
    with open(statsfile, 'r') as inFile:
        for line in inFile:
            if line.startswith('P10 ') or line.startswith('P_10 '):
                tokens = line.strip().split()
                p10s.append((tokens[1], tokens[2]))
            if line.startswith('ndcg_cut_1000 '):
                tokens = line.strip().split()
                ndcgs.append((tokens[1], tokens[2]))
            if line.startswith('num_rel_ret '):
                tokens = line.strip().split()
                num_rel_rets.append((tokens[1], tokens[2]))
            if line.startswith('map_cut_1000 '):
                tokens = line.strip().split()
                maps.append((tokens[1], tokens[2]))
            if line.startswith('recall_1000 '):
                tokens = line.strip().split()
                recalls.append((tokens[1], tokens[2]))
            if line.startswith('num_rel '):
                tokens = line.strip().split()
                num_rels.append((tokens[1], tokens[2]))
            if line.startswith('Rprec '):
                tokens = line.strip().split()
                rprecs.append((tokens[1], tokens[2]))

    with open(statsfile2, 'r') as inFile2:
        for line in inFile2:
            if line.startswith('infNDCG'):
                tokens = line.strip().split()
                infndcgs.append((tokens[1], tokens[2]))

        p10s.sort(key=key)
        num_rel_rets.sort(key=key)
        num_rels.sort(key=key)
        ndcgs.sort(key=key)
        maps.sort(key=key)
        recalls.sort(key=key)
        rprecs.sort(key=key)
        infndcgs.sort(key=key)

    re5000s = []
    with open(qfilename, 'a+') as qstatsfile:
        qstatsfile.write(os.path.basename(filenamebase))
        for p10, recall, num_rel, num_rel_ret, ndcg, mapp, rprec, infNDCG in zip(p10s, recalls, num_rels, num_rel_rets, ndcgs, maps, rprecs, infndcgs):
            qno = p10[0]
            if qno != 'all':
                re5000 = int(num_rel_ret[1])/int(num_rel[1])
                re5000s.append(re5000)
            else:
                re5000 = sum(re5000s)/len(re5000s)
            qstatsfile.write(',{},{},{},{},{},{},{}'.format(p10[1], recall[1], re5000, ndcg[1], mapp[1], rprec[1], infNDCG[1]))

        qstatsfile.write('\n')


# Create average stats file for a given TREC format 
def main(filename, outfilename='p10_stats.txt', avgstatsfile='avgstats.csv'):
    filenamebase = filename[:-8]
    stats = ('P_10', 'recall_1000', 'ndcg', 'map', 'num_ret', 'num_rel_ret', 'iprec_at_recall_0.50', 'ndcg_cut_1000',
             'map_cut_1000')

    if not os.path.exists(avgstatsfile):
        with open(avgstatsfile, 'w') as avgoutfile:
            avgoutfile.write('formulation,' + ','.join(stats) + '\n')

    statsfile, statsfile2 = create_stats(filename)

    p10s = []
    ndcgs = []
    reprecs = []
    avgstats = {}
    with open(statsfile, 'r') as inFile:
        with open(outfilename, 'a+') as p10File, open('NDCG_stats.txt', 'a+') as ndcgFile, open('rprec_stats.txt', 'a+') as rprecFile:
            p10File.write(filenamebase + '\n')
            ndcgFile.write(filenamebase + '\n')
            rprecFile.write(filenamebase + '\n')
            for line in inFile:
                if line.startswith('P10 ') or line.startswith('P_10 '):
                    tokens = line.strip().split()
                    p10s.append((tokens[1], tokens[2]))
                if line.startswith('ndcg '):
                    tokens = line.strip().split()
                    ndcgs.append((tokens[1], tokens[2]))
                if line.startswith('Rprec '):
                    tokens = line.strip().split()
                    reprecs.append((tokens[1], tokens[2]))

                if line.startswith(stats):
                    tokens = line.strip().split()
                    if tokens[1] == 'all':
                        avgstats[tokens[0]] = tokens[2]

            p10s.sort(key=key)
            ndcgs.sort(key=key)
            reprecs.sort(key=key)

            for p10, ndcg, map in zip(p10s, ndcgs, reprecs):
                p10File.write(p10[1] + '\n')
                ndcgFile.write(ndcg[1] + '\n')
                rprecFile.write(map[1] + '\n')

    with open(avgstatsfile, 'a+') as outfile:
        outfile.write(filenamebase + ',' + ','.join([avgstats[stat] for stat in stats]) + '\n')


# sort keys
def key(item):
    qno = item[0]
    if qno.isdigit():
        return int(qno)
    else:
        return 10000

def strkey(item):
    qno = item
    if qno.isdigit():
        return int(qno)
    else:
        return 10000

if __name__ == '__main__':
    main(sys.argv[1])
