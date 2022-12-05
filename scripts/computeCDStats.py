import argparse
import dendropy
import os
from pathlib import Path
import glob
import pickle

'''
Below function is from the Quintet Rooting repository at
https://github.com/ytabatabaee/Quintet-Rooting/blob/main/scripts/clade_distance.py
written by Yasamin Tabatabaee, based on code from NJMerge software
written by Erin K. Molloy, released under a 3-clause BSD license
(see https://opensource.org/licenses/BSD-3-Clause).
'''
def clade_distance(tr1, tr2):
    from dendropy.calculate.treecompare \
        import symmetric_difference

    lb1 = set([l.taxon.label for l in tr1.leaf_nodes()])
    lb2 = set([l.taxon.label for l in tr2.leaf_nodes()])

    com = lb1.intersection(lb2)
    if com != lb1 or com != lb2:
        com = list(com)
        tns = dendropy.TaxonNamespace(com)

        tr1.retain_taxa_with_labels(com)
        tr1.migrate_taxon_namespace(tns)

        tr2.retain_taxa_with_labels(com)
        tr2.migrate_taxon_namespace(tns)
    com = list(com)

    tr1.update_bipartitions()
    tr2.update_bipartitions()

    nl = len(com)
    cd = symmetric_difference(tr1, tr2) / (2*nl - 4)

    return nl, cd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compute clade distances of data')

    parser.add_argument('-c', '--computedtrees', type=str, required=True,
                        help='Computed trees folder')
    parser.add_argument('-n', '--name', type=str, required=True,
                        help='Name of computed tree files')
    parser.add_argument('-t', '--truetrees', type=str, required=True,
                        help='True (rooted) species trees folder')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='File to write output data')
    parser.add_argument('--taxa', type=str, required=True, help='List of taxa datasets')
    
    args = parser.parse_args()

    datasets = list([str(i)+'_taxon' for i in args.taxa.split(',')])
    gteeList = ['true', '1500', '1000', '500', '250']
    nReplicates = 20
    try:
        with open(args.output, 'rb') as outputFile:
            outputData = pickle.load(outputFile)
    except:
        outputData = {}
    for datasetNum,dataset in enumerate(datasets):
        trueTreePath = os.path.join(args.truetrees, dataset)
        truePathList = list(Path().glob(trueTreePath + '/**/model-species.tre'))

        outputData[dataset] = {}
        for gtee in gteeList:
            outputData[dataset][gtee] = [0.]*nReplicates
            for rep in range(1, nReplicates+1):
                for truePath in truePathList:
                    iterateNum = truePath.parent.stem
                    computedPath = os.path.join(args.computedtrees,
                                                 dataset, iterateNum, gtee,
                                                 'R'+str(rep), args.name)

                    taxa = dendropy.TaxonNamespace()
                    trueTree = dendropy.Tree.get(path=str(truePath),
                                                schema='newick',
                                                taxon_namespace=taxa,
                                                rooting='force-rooted')
                    computedTree = dendropy.Tree.get(path=str(computedPath),
                                                schema='newick',
                                                taxon_namespace=taxa,
                                                rooting='force-rooted')
                    _, cd = clade_distance(trueTree, computedTree)
                    outputData[dataset][gtee][rep-1] += cd/len(truePathList)

        print(f'\rDone with {datasetNum+1}/{len(datasets)}..', end='')
    print()
    with open(args.output, 'wb') as outputFile:
        pickle.dump(outputData, outputFile)

