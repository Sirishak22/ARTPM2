######################################
# Quick and dirty grid searches for LeToR hyperparameters, not intended for extensive use
# Lowell Milliken
######################################
import validation as val
import extract_stat as es


def do_search():
    # ranker = 'Random Forests'
    ranker = 'ListNet'
    stats = ('P_10', 'recall_1000', 'ndcg_cut_1000', 'map_cut_1000')
    outfilename = 'LN_grid_search_extra.csv'

    with open(outfilename, 'w') as outfile:
        outfile.write('learning rate,epoch,norms,' + ','.join(stats) + '\n')
    #     outfile.write('frate,tree,leaf,shrinkage,norm,' + ','.join(stats) + '\n')

    # frate = 0.2
    nscores = False
    for norm in ('linear', 'none'):
        # for frate in (0.3, 0.5, 0.8):
        #     for tree in (1,5,10):
        #         for leaf in (50, 100, 150, 200):

        # tree = 1
        # leaf = 50
        for lr in ('0.000001', '0.00001', '0.0001', '0.001', '0.01', '0.1', '1', '10', '100'):
            for epoch in (1500, 3000, 4500):
                    # params = {'-frate': frate, '-tree': tree, '-leaf': leaf, '-shrinkage': '0.01'}
                    params = {'-lr': lr, '-epoch': epoch}
                    if norm != 'none':
                        params['-norm'] = norm
                    runfilename = val.do_cv(meta=False, splitdrugs=False, textlen=True, indriscore=False,
                                            ranker=ranker, rparams=params, normscores=nscores)

                    statsfile = es.create_stats(runfilename)
                    avgstats = es.get_stats(statsfile, stats)

                    with open(outfilename, 'a+') as outfile:
                        # outfile.write('{},{},{},{},{},'.format(frate, tree, leaf, 0.01, norm) + ','.join([avgstats[stat] for stat in stats]) + '\n')
                        outfile.write('{},{},{},'.format(lr, epoch, norm) + ','.join(
                            [avgstats[stat] for stat in stats]) + '\n')


if __name__ == '__main__':
    do_search()