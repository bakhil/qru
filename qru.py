'''
Most of this code is from Quintet Rooting developed by
Yasamin Tabatabaee (see https://github.com/ytabatabaee/Quintet-Rooting).
Quintet Rooting was released under a GNU GPL v3.0 license.
'''
import argparse
import time
import dendropy
import numpy as np
import sys
from table_five import TreeSet
import subprocess
import re

sys.path.append('Quintet_Rooting')
from qr.adr_theory import *
from qr.fitness_cost import *
from qr.quintet_sampling import *
from qr.utils import *

# !! Define name of root which will be inserted
NAME_OF_ROOT = 'ROOTLEAFTOINSERT'

def main(args):
    st_time = time.time()
    # !! Change path to script for reading ADR quintets
    script_path = os.path.join(os.path.realpath(__file__).rsplit("/", 1)[0], 'Quintet_Rooting')

    # input args
    species_tree_path = args.speciestree
    gene_tree_path = args.genetrees
    output_path = args.outputtree
    sampling_method = args.samplingmethod.lower()
    random.seed(args.seed)
    cost_func = args.cost.lower()
    shape_coef = args.coef
    mult_le = args.multiplicity
    abratio = args.abratio

    ## !! folder for writing junk files - and other filenames
    #junkFolder = args.junk
    #astral_input_filename = os.path.join(junkFolder, 'astral_input.tre')
    #astral_bipartition_filename = os.path.join(junkFolder, 'astral_bipartitions.tre')

    # !! Define the path to ASTRAL executable
    astral_path = args.astralpath

    # !! Commented out header
    '''
    header = """*********************************
*     Quintet Rooting """ + __version__ + """    *
*********************************"""
    sys.stdout.write(header + '\n')
    '''


    tns = dendropy.TaxonNamespace()
    unrooted_species = dendropy.Tree.get(path=species_tree_path, schema='newick',
                                         taxon_namespace=tns, rooting="force-unrooted", suppress_edge_lengths=True)
    if len(tns) < 5:
        raise Exception("Species tree " + species_tree_path + " has less than 5 taxa!\n")
    gene_trees = TreeSet(gene_tree_path)

    # reading fixed quintet topology files
    tns_base = dendropy.TaxonNamespace()
    unrooted_quintets_base = dendropy.TreeList.get(path=script_path + '/qr/topologies/quintets.tre',
                                                   taxon_namespace=tns_base, schema='newick')
    rooted_quintets_base = dendropy.TreeList(taxon_namespace=tns_base)
    rooted_quintets_base.read(path=script_path + '/qr/topologies/caterpillar.tre', schema='newick',
                              rooting="default-rooted")
    rooted_quintets_base.read(path=script_path + '/qr/topologies/pseudo_caterpillar.tre', schema='newick',
                              rooting="default-rooted")
    rooted_quintets_base.read(path=script_path + '/qr/topologies/balanced.tre', schema='newick',
                              rooting="default-rooted")
    rooted_quintet_indices = np.load(script_path + '/qr/rooted_quintet_indices.npy')

    #sys.stdout.write('Loading time: %.2f sec\n' % (time.time() - st_time))
    ss_time = time.time()

    # !! Comment this out here and add for each sampled quintet
    ## search space of rooted trees
    #rooted_candidates = get_all_rooted_trees(unrooted_species)
    #r_score = np.zeros(len(rooted_candidates))

    #sys.stdout.write('Creating search space time: %.2f sec\n' % (time.time() - ss_time))
    sm_time = time.time()

    # set of sampled quintets
    taxon_set = [t.label for t in tns]
    sample_quintet_taxa = []
    if len(taxon_set) == 5 or sampling_method == 'exh':
        sample_quintet_taxa = list(itertools.combinations(taxon_set, 5))
    elif sampling_method == 'tc':
        sample_quintet_taxa = triplet_cover_sample(taxon_set)
    elif sampling_method == 'le':
        sample_quintet_taxa = linear_quintet_encoding_sample(unrooted_species, taxon_set, mult_le)
    elif sampling_method == 'rl':
        sample_quintet_taxa = random_linear_sample(taxon_set)

    #sys.stdout.write('Quintet sampling time: %.2f sec\n' % (time.time() - sm_time))
    proc_time = time.time()

    # !! Comment out output
    #sys.stdout.write("Number of taxa (n): %d\n" % len(tns))
    #sys.stdout.write("Number of gene trees (k): %d\n" % len(gene_trees))
    #sys.stdout.write("Size of search space (|R|): %d\n" % len(rooted_candidates))
    #sys.stdout.write("Size of sampled quintets set (|Q*|): %d\n" % len(sample_quintet_taxa))

    # preprocessing
    quintet_scores = np.zeros((len(sample_quintet_taxa), 7))
    quintet_unrooted_indices = np.zeros(len(sample_quintet_taxa), dtype=int)
    quintets_r_all = []
    rooted_quintets_from_samples = []

    for j in range(len(sample_quintet_taxa)):
        q_taxa = sample_quintet_taxa[j]
        # !! Find rooted quintets - below three lines added
        # !! (first line is just defintion of subtree_u from previous code)
        # search space of rooted trees
        subtree_u = unrooted_species.extract_tree_with_taxa_labels(labels=q_taxa, suppress_unifurcations=True)
        rooted_candidates = get_all_rooted_trees(subtree_u)
        r_score = np.zeros(len(rooted_candidates))

        quintets_u = [
            dendropy.Tree.get(data=map_taxon_namespace(str(q), q_taxa) + ';', schema='newick', rooting='force-unrooted',
                              taxon_namespace=tns) for q in unrooted_quintets_base]
        quintets_r = [
            dendropy.Tree.get(data=map_taxon_namespace(str(q), q_taxa) + ';', schema='newick', rooting='force-rooted',
                              taxon_namespace=tns) for q in rooted_quintets_base]
        #subtree_u = unrooted_species.extract_tree_with_taxa_labels(labels=q_taxa, suppress_unifurcations=True)
        quintet_counts = np.asarray(gene_trees.tally_single_quintet(q_taxa))
        quintet_normalizer = sum(quintet_counts) if args.normalized else len(gene_trees)
        quintet_tree_dist = quintet_counts
        if quintet_normalizer != 0:
            quintet_tree_dist = quintet_tree_dist / quintet_normalizer
        quintet_unrooted_indices[j] = get_quintet_unrooted_index(subtree_u, quintets_u)
        quintet_scores[j] = compute_cost_rooted_quintets(quintet_tree_dist, quintet_unrooted_indices[j],
                                                         rooted_quintet_indices, cost_func, len(gene_trees),
                                                         len(sample_quintet_taxa), shape_coef, abratio)

        # !! Comment out below line since it's not required
        #quintets_r_all.append(quintets_r)

        # !! Find the best quintet here
        # Code from below - see commented lines below for changes
        # computing scores
        min_score = sys.maxsize
        for i in range(len(rooted_candidates)):
            r = rooted_candidates[i]
            #for j in range(len(sample_quintet_taxa)):
            #q_taxa = sample_quintet_taxa[j]
            subtree_r = r.extract_tree_with_taxa_labels(labels=q_taxa, suppress_unifurcations=True)
            r_idx = get_quintet_rooted_index(subtree_r, quintets_r, quintet_unrooted_indices[j])
            r_score[i] += quintet_scores[j][r_idx]
            #if not args.confidencescore and r_score[i] > min_score:
            #    break
            if r_score[i] < min_score:
                min_score = r_score[i]

        min_idx = np.argmin(r_score)
        # !! Comment out below lines and write rooted quintet for ASTRAL
        #with open(output_path, 'w') as fp:
        #    fp.write(str(rooted_candidates[min_idx]) + ';\n')
        
        # !! Store rooted quintet unrooted using new root
        rooted_quintets_from_samples.append(rooted_candidates[min_idx])

    # !! Write to input and bipartition files for ASTRAL
    astral_input_str = ''
    for tree in rooted_quintets_from_samples:
        outputStr = '(' + str(tree) + ',' + NAME_OF_ROOT + ');\n'
        astral_input_str += outputStr
    astral_input_str += (str(unrooted_species) + ';\n')
    #with open(astral_input_filename, 'w') as astral_input_file:
    #    for tree in rooted_quintets_from_samples:
    #        outputStr = '(' + str(tree) + ',' + NAME_OF_ROOT + ');\n'
    #        astral_input_file.write(outputStr)
    #    astral_input_file.write(str(unrooted_species) + ';\n')
    astral_bipartition_str = ''
    for taxon_label in taxon_set:
        new_tree = unrooted_species.clone(depth=1)
        root_at_node = new_tree.find_node_with_taxon_label(taxon_label)
        new_tree.reroot_at_edge(root_at_node.edge, update_bipartitions=True)
        outputStr = '(' + str(new_tree) + ',' + NAME_OF_ROOT + ');\n'
        astral_bipartition_str += outputStr
    #with open(astral_bipartition_filename, 'w') as astral_bipartition_file:
    #    for taxon_label in taxon_set:
    #        new_tree = unrooted_species.clone(depth=1)
    #        root_at_node = new_tree.find_node_with_taxon_label(taxon_label)
    #        new_tree.reroot_at_edge(root_at_node.edge, update_bipartitions=True)
    #        outputStr = '(' + str(new_tree) + ',' + NAME_OF_ROOT + ');\n'
    #        astral_bipartition_file.write(outputStr)
    
    # !! Run ASTRAL to get output
    cmd = 'java -jar ' + astral_path + ' -i <(printf "' + astral_input_str + '") ' + \
                    '-f <(printf "' + astral_bipartition_str + '")'
    cmd = 'bash -c \'' + cmd + '\''
    ps = subprocess.run(cmd, shell=True, capture_output=True)
    output = ps.stderr.decode(sys.stderr.encoding)
    matches = re.findall(r'\n(\(.*\);)\n', output)
    if len(matches) == 0:
        raise Exception('ASTRAL failed to run properly!!')

    # !! Get tree and write to file
    outputTreeStr = matches[0]
    newtxn = dendropy.TaxonNamespace()
    outputTree = dendropy.Tree.get(data=outputTreeStr, schema='newick',
                                rooting='force-rooted', taxon_namespace=newtxn)
    nodeToRoot = outputTree.find_node_with_taxon_label(NAME_OF_ROOT)
    outputTree.reroot_at_edge(nodeToRoot.edge, update_bipartitions=True)
    finalOutput = outputTree.extract_tree_with_taxa_labels(taxon_set)
    with open(output_path, 'w') as fp:
        fp.write(str(finalOutput) + ';\n')

    #sys.stdout.write('Preprocessing time: %.2f sec\n' % (time.time() - proc_time))
    sc_time = time.time()

    # !! Comment out below lines here and add to each quintet
    ## computing scores
    #min_score = sys.maxsize
    #for i in range(len(rooted_candidates)):
    #    r = rooted_candidates[i]
    #    for j in range(len(sample_quintet_taxa)):
    #        q_taxa = sample_quintet_taxa[j]
    #        subtree_r = r.extract_tree_with_taxa_labels(labels=q_taxa, suppress_unifurcations=True)
    #        r_idx = get_quintet_rooted_index(subtree_r, quintets_r_all[j], quintet_unrooted_indices[j])
    #        r_score[i] += quintet_scores[j][r_idx]
    #        if not args.confidencescore and r_score[i] > min_score:
    #            break
    #    if r_score[i] < min_score:
    #        min_score = r_score[i]

    #min_idx = np.argmin(r_score)
    #with open(output_path, 'w') as fp:
    #    fp.write(str(rooted_candidates[min_idx]) + ';\n')

    #sys.stdout.write('Scoring time: %.2f sec\n' % (time.time() - sc_time))
    #sys.stdout.write('Best rooting: \n%s \n' % str(rooted_candidates[min_idx]))

    ## computing confidence scores
    #if args.confidencescore:
    #    sys.stdout.write('Scores of all rooted trees:\n %s \n' % str(r_score))
    #    confidence_scores = (np.max(r_score) - r_score) / np.sum(np.max(r_score) - r_score)
    #    tree_ranking_indices = np.argsort(r_score)
    #    with open(output_path + ".rank.cfn", 'w') as fp:
    #        for i in tree_ranking_indices:
    #            fp.write(str(rooted_candidates[i]) + ';\n')
    #            fp.write(str(confidence_scores[i]) + '\n')

    #sys.stdout.write('Total execution time: %.2f sec\n' % (time.time() - st_time))


def compute_cost_rooted_quintets(u_distribution, u_idx, rooted_quintet_indices, cost_func, k, q_size, shape_coef, abratio):
    """
    Scores the 7 possible rootings of an unrooted quintet
    :param np.ndarray u_distribution: unrooted quintet tree probability distribution
    :param int u_idx: index of unrooted binary tree
    :param np.ndarray rooted_quintet_indices: indices of partial orders for all rooted quintet trees
    :param str cost_func: type of the fitness function
    :rtype: np.ndarray
    """
    rooted_tree_indices = u2r_mapping[u_idx]
    costs = np.zeros(7)
    for i in range(7):
        idx = rooted_tree_indices[i]
        unlabeled_topology = idx_2_unlabeled_topology(idx)
        indices = rooted_quintet_indices[idx]
        costs[i] = cost(u_distribution, indices, unlabeled_topology, cost_func, k, q_size, shape_coef, abratio)
    return costs


def get_all_rooted_trees(unrooted_tree):
    """
    Generates all the possible rooted trees with a given unrooted topology
    :param dendropy.Tree unrooted_tree: an unrooted tree topology
    :rtype: list
    """
    rooted_candidates = []
    tree = dendropy.Tree(unrooted_tree)
    for edge in tree.preorder_edge_iter():
        try:
            tree.reroot_at_edge(edge, update_bipartitions=True)
            rooted_candidates.append(dendropy.Tree(tree))
        except:
            continue
    # removing duplicates
    rooted_candidates[0].resolve_polytomies(update_bipartitions=True)
    for i in range(1, len(rooted_candidates)):
        if dendropy.calculate.treecompare.symmetric_difference(rooted_candidates[0], rooted_candidates[i]) == 0:
            rooted_candidates.pop(0)
            break
    return rooted_candidates


def parse_args():
    # !! Change header
    parser = argparse.ArgumentParser(description='Quintet Rooting via unrooting')

    parser.add_argument("-t", "--speciestree", type=str,
                        help="input unrooted species tree in newick format",
                        required=True, default=None)

    parser.add_argument("-g", "--genetrees", type=str,
                        help="input gene trees in newick format",
                        required=True, default=None)

    parser.add_argument("-o", "--outputtree", type=str,
                        help="output file containing a rooted species tree",
                        required=True, default=None)

    parser.add_argument("-sm", "--samplingmethod", type=str,
                        help="quintet sampling method (LE for linear encoding (default), EXH for exhaustive",
                        required=False, default='LE')

    parser.add_argument("-c", "--cost", type=str,
                        help="cost function (STAR for running QR-STAR)",
                        required=False, default='d')

    parser.add_argument("-cfs", "--confidencescore", action='store_true',
                        help="output confidence scores for each possible rooted tree as well as a ranking")

    parser.add_argument("-mult", "--multiplicity", type=int,
                        help="multiplicity (number of quintets mapped to each edge) in QR-LE",
                        required=False, default=1)

    parser.add_argument("-norm", "--normalized", action='store_true',
                        help="normalization for unresolved gene trees or missing taxa",
                        required=False, default=False)

    parser.add_argument("-coef", "--coef", type=float,
                        help="coefficient for shape penalty term in QR-STAR", required=False, default=0)

    parser.add_argument("-abratio", "--abratio", type=float,
                        help="Ratio between invariant and inequality penalties used in QR-STAR", required=False, default=1)

    parser.add_argument("-rs", "--seed", type=int,
                        help="random seed", required=False, default=1234)

    parser.add_argument('-j', '--junk', type=str,
                        help='Folder for writing temporary files',
                        required=False, default='junk')

    parser.add_argument('-apath', '--astralpath', type=str,
                        default='ASTRAL-modified/astral.5.7.3.modified.jar', required=False)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main(parse_args())
