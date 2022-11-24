import argparse
import dendropy
import os
from pathlib import Path

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
    parser = argparse.ArgumentParser(description='Find the best rooting of estimated species trees.')
    parser.add_argument('-s', '--speciestree', type=str, required=True,
                        help='Folder containing estimated unrooted species trees')
    parser.add_argument('-t', '--truetrees', type=str, required=True,
                        help='True (rooted) species trees')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Folder to write output trees')
    
    args = parser.parse_args()

    for folder, _, fileNames in os.walk(args.truetrees):
        for fileName in fileNames:
            _, ext = os.path.splitext(fileName)
            if ext != '.tre' and ext != '.tree':
                continue
            trueTreePath = Path(os.path.join(folder, fileName))
            taxa = dendropy.TaxonNamespace()
            trueTree = dendropy.Tree.get(path=str(trueTreePath),
                                         schema='newick',
                                         taxon_namespace=taxa,
                                         rooting='force-rooted')

            nPaths = len(Path(args.truetrees).parts)
            speciesTreesFolder = os.path.join(args.speciestree,
                                str(Path(*trueTreePath.parts[nPaths:-1])))
            outputTreesFolder = os.path.join(args.output,
                                str(Path(*trueTreePath.parts[nPaths:-1])))
            nPaths = len(Path(speciesTreesFolder).parts)
            for folder, _, speciesFileNames in os.walk(speciesTreesFolder):
                for speciesFileName in speciesFileNames:
                    speciesTreePath = Path(folder) / speciesFileName
                    outputTreePath = Path(outputTreesFolder) / \
                                Path(*Path(speciesTreePath).parts[nPaths:-1]) / \
                                ('best-root-' + speciesFileName)
            
                    speciesTree = dendropy.Tree.get(path=str(speciesTreePath),
                                            schema='newick',
                                            taxon_namespace=taxa,
                                            rooting='force-unrooted')
            
            
                    bestCd = 1e10
                    bestTree = speciesTree.clone(depth=1)
                    numEdges = 2*len(taxa) - 3
                    for edge in list(speciesTree.preorder_edge_iter())[-numEdges:]:
                        speciesTree.reroot_at_edge(edge, update_bipartitions=True)
                        nl, cd = clade_distance(speciesTree, trueTree)
                        if cd < bestCd:
                            bestCd = cd
                            bestTree = speciesTree.clone(depth=1)
                        speciesTree.deroot()

                    outputTreePath.parent.mkdir(parents=True, exist_ok=True)
                    bestTree.write_to_path(str(outputTreePath), schema='newick')

