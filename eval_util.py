##############################################
# Methods to help evaluate ARtPM performance.
# Lowell Milliken
##############################################

import validation as val
import extract_stat as es
from scipy import stats
from itertools import combinations


# Do training and validation on different Retrieval stage formulations.
def do_val(indri):
    if not indri:
        file = val.do_cv(
            unknown_docs_filename='topics2017_m_as_ex_tr_fd_ft_nsh_d_ds_lnm_noexp_prf-2-20-0.5-0.5_large_gfix_alldocs.txt',
            rparams={'-lr': 0.1, '-epoch': 3000}, scorefile=None)
        es.main(file)
        file = val.do_cv(unknown_docs_filename='topics2017_m_as_ex_tr_fd_ft_nsh_d_ds_lnm_large_gfix_alldocs.txt',
                   rparams={'-lr': 0.1, '-epoch': 3000}, scorefile=None)
        es.main(file)
        file = val.do_cv(unknown_docs_filename='topics2017_m_as_ex_tr_fd_ft_nsh_lnm_prf-2-20-0.5-0.5_large_gfix_alldocs.txt',
                  rparams={'-lr': 0.1, '-epoch': 3000}, scorefile=None)
        es.main(file)
        file = val.do_cv(unknown_docs_filename='topics2017_m_as_ex_tr_fd_ft_nsh_d_ds_prf-2-20-0.5-0.5_large_gfix_alldocs.txt',
                  rparams={'-lr': 0.1, '-epoch': 3000}, scorefile=None)
        es.main(file)
        file = val.do_cv(unknown_docs_filename='topics2017_m_as_ex_fd_ft_nsh_d_ds_lnm_prf-2-20-0.5-0.5_large_gfix_alldocs.txt',
                  rparams={'-lr': 0.1, '-epoch': 3000}, scorefile=None)
        es.main(file)
        file = val.do_cv(unknown_docs_filename='topics2017_m_as_tr_fd_ft_nsh_d_ds_lnm_prf-2-20-0.5-0.5_large_gfix_alldocs.txt',
                  rparams={'-lr': 0.1, '-epoch': 3000}, scorefile=None)
        es.main(file)
        file = val.do_cv(unknown_docs_filename='topics2017_m_as_ex_tr_nd_ft_nsh_d_ds_lnm_prf-2-20-0.5-0.5_large_gfix_alldocs.txt',
                  rparams={'-lr': 0.1, '-epoch': 3000}, scorefile=None)
        es.main(file)
        file = val.do_cv(unknown_docs_filename='topics2017_m_as_ex_tr_fd_ft_nsh_d_ds_lnm_prf-2-20-0.5-0.5_large_gfix_alldocs.txt',
                  rparams={'-lr': 0.1, '-epoch': 3000}, scorefile=None)
        es.main(file)
    else:
        file = val.do_cv(unknown_docs_filename='topics2017_m_as_ex_tr_fd_ft_nsh_d_ds_lnm_large_gfix_alldocs.txt',
                         rparams={'-lr': 0.1, '-epoch': 3000}, indriscore=True, scorefile=None)
        es.main(file)
        file = val.do_cv(unknown_docs_filename='topics2017_m_as_ex_tr_fd_ft_nsh_d_ds_lnm_noexp_prf-2-20-0.5-0.5_large_gfix_alldocs.txt',
                         rparams={'-lr': 0.1, '-epoch': 3000}, indriscore=True, scorefile=None)
        es.main(file)
        file = val.do_cv(unknown_docs_filename='topics2017_m_as_ex_tr_fd_ft_nsh_lnm_prf-2-20-0.5-0.5_large_gfix_alldocs.txt',
                  rparams={'-lr': 0.1, '-epoch': 3000}, indriscore=True, scorefile=None)
        es.main(file)
        file = val.do_cv(unknown_docs_filename='topics2017_m_as_ex_tr_fd_ft_nsh_d_ds_prf-2-20-0.5-0.5_large_gfix_alldocs.txt',
                  rparams={'-lr': 0.1, '-epoch': 3000}, indriscore=True, scorefile=None)
        es.main(file)
        file = val.do_cv(unknown_docs_filename='topics2017_m_as_ex_fd_ft_nsh_d_ds_lnm_prf-2-20-0.5-0.5_large_gfix_alldocs.txt',
                  rparams={'-lr': 0.1, '-epoch': 3000}, indriscore=True, scorefile=None)
        es.main(file)
        file = val.do_cv(unknown_docs_filename='topics2017_m_as_tr_fd_ft_nsh_d_ds_lnm_prf-2-20-0.5-0.5_large_gfix_alldocs.txt',
                  rparams={'-lr': 0.1, '-epoch': 3000}, indriscore=True, scorefile=None)
        es.main(file)
        file = val.do_cv(unknown_docs_filename='topics2017_m_as_ex_tr_nd_ft_nsh_d_ds_lnm_prf-2-20-0.5-0.5_large_gfix_alldocs.txt',
                  rparams={'-lr': 0.1, '-epoch': 3000}, indriscore=True, scorefile=None)
        es.main(file)
        file = val.do_cv(unknown_docs_filename='topics2017_m_as_ex_tr_fd_ft_nsh_d_ds_lnm_prf-2-20-0.5-0.5_large_gfix_alldocs.txt',
                  rparams={'-lr': 0.1, '-epoch': 3000}, indriscore=True, scorefile=None)
        es.main(file)


# Run t-tests on all combinations of runs in "runs" for "stat" metric.
# "all_stats" contains query level metrics for all runs.
def do_t(runs, all_stats, stat):
    results = {}

    for pair in combinations(runs, 2):
        results[pair] = stats.ttest_rel(all_stats[pair[0]][stat], all_stats[pair[1]][stat])[1]

    return results


# save all stats in check_stats to outfile csv
def run_t(runs, outfile, check_stats, all_stats):
    results = {}
    for stat in check_stats:
        results[stat] = do_t(runs, all_stats, stat)
        outfile.write('\n{}\n'.format(stat))
        outfile.write(',' + ','.join(runs))
        for run1 in runs:
            outfile.write('\n{},'.format(run1))
            for run2 in runs:
                tup = (run1, run2)
                if tup in results[stat]:
                    outfile.write('{},'.format(results[stat][tup]))
                else:
                    outfile.write('---,')

    return results

# do t-tests for the listed runs and the listed stats.
def do_all_t(outfilename='ttestpscores4.csv'):
    check_stats = ('P_10', 'P_30', 'P_100', 'recall_1000', 'num_rel_ret', 'num_rel', 'ndcg_cut_1000', 'map_cut_1000',
                   'Rprec_mult_0.20', 'Rprec_mult_0.40', 'Rprec_mult_0.60', 'Rprec_mult_0.80', 'Rprec_mult_1.00',
                   'Rprec_mult_1.20', 'Rprec_mult_1.40', 'Rprec_mult_1.60', 'Rprec_mult_1.80', 'Rprec_mult_2.00',
                   'infNDCG', 'iprec_at_recall_0.00', 'iprec_at_recall_0.10', 'iprec_at_recall_0.20',
                   'iprec_at_recall_0.30', 'iprec_at_recall_0.40', 'iprec_at_recall_0.50', 'iprec_at_recall_0.60',
                   'iprec_at_recall_0.70', 'iprec_at_recall_0.80', 'iprec_at_recall_0.90', 'iprec_at_recall_1.00',
                   'Recall_5000')
    all_stats = load_stats(check_stats)

    with open(outfilename, 'w') as outfile:
        runs1 = ['PubMed search', 'Indri baseline']

        runs2 = ['no letor no term expansion', 'no letor no prf', 'no letor no drugs',
                 'no letor no mutations', 'no letor no treatment words', 'no letor no exact phrases',
                 'no letor no demographic', 'no letor no ablation']

        runs3 = ['no term expansion', 'no prf', 'no drugs', 'no mutations', 'no treatment words', 'no exact phrases',
                 'no demographic', 'no ablation']

        runs4 = ['letor with indri no prf', 'letor with indri no term expansion', 'letor with indri no drugs',
                 'letor with indri no mutations', 'letor with indri no treatment words', 'letor with indri no exact phrases',
                 'letor with indri no demographic', 'letor with indri no ablation']

        runs5 = ['letor with high prec scores no prf', 'letor with high prec scores no treatment words',
                 'letor with high prec scores no drugs', 'letor with high prec scores no mutations',
                 'letor with high prec scores no term expansion', 'letor with high prec scores no exact phrases',
                 'letor with high prec scores no demographic', 'letor with high prec scores no ablation']

        runs6 = ['ARtPM without query formulation', 'ARtPM without external drug knowledge results', 'ARtPM without LETOR', 'ARtPM system']

        runs7 = ['mugpubbase', 'mugpubboost', 'UTDHLTAF', 'UTDHLTFF', 'UD_GU_SA_5']

        runs = runs1 + runs2 + runs3 + runs4 + runs5 + runs6 + runs7
        runs = ['ARtPM without LETOR', 'no exact phrases']
        results = run_t(runs, outfile, check_stats, all_stats)

    return results


_basestats = ['P_10', 'Recall', 'num_rel_ret', 'ndcg', 'map', 'Rprec', 'infNDCG']


# convert csv to latex format for tables
def to_latex(check_stats, filename='querystats5.csv', outfilename='avgstats2.tx'):
    all_stats = load_stats(check_stats, filename, avg=True)
    with open(outfilename, 'w') as outfile:
        outfile.write('formulation & ' + '&'.join(check_stats) + '\\\\\n')
        for form, fstats in all_stats.items():
            outfile.write(form)
            for stat in check_stats:
                outfile.write('&{:0.3f}'.format(fstats[stat]))

            outfile.write('\\\\\n')


# load query stats from a csv stat file.
def load_stats(check_stats, statsfilename='querystats6.csv', avg=False, querycount=30):
    fstats = {}
    with open(statsfilename, 'r') as statsfile:
        statsfile.readline()
        for line in statsfile:
            tokens = line.split(',')
            fstats[tokens[0]] = dict()
            for stat in check_stats:
                fstats[tokens[0]][stat] = []
            position = 1
            for i in range(querycount):
                for stat in check_stats:
                    fstats[tokens[0]][stat].append(float(tokens[position]))
                    position += 1

            if avg:
                for stat in check_stats:
                    fstats[tokens[0]][stat] = float(tokens[position])
                    position += 1

    return fstats
