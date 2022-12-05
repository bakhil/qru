import argparse
import os
from pathlib import Path

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run QRU over QR-paper datasets')
    parser.add_argument('-t', '--speciestree', type=str, required=True,
                        help='Folder containing estimated unrooted species trees')
    parser.add_argument('-sn', '--speciesname', type=str, required=True,
                        help='Name of species files')
    parser.add_argument('-g', '--genetrees', type=str, required=True,
                        help='Folder containing (estimated) gene trees')
    parser.add_argument('-gn', '--genename', type=str, required=True,
                        help='Name of gene files')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Folder to write output trees')
    parser.add_argument('-on', '--outputname', type=str, required=True,
                        help='Name of output files')
    parser.add_argument('--taxa', type=str, required=True,
                        help='Comma separated numbers of taxa to use')
    parser.add_argument('--gtee', type=str, required=True,
                        help='Comma separated list of GTEE conditions')
    parser.add_argument('--qpath', type=str, required=False,
                        default='../qru.py', help='Path to qru.py')
    
    args = parser.parse_args()
    datasets = [t+'_taxon' for t in args.taxa.split(',')]
    gteeList = args.gtee.split(',')
    nReplicates = 20
    for dataset in datasets:
        speciesFolder = os.path.join(args.speciestree, dataset)
        geneFolder = os.path.join(args.genetrees, dataset)
        outputFolder = os.path.join(args.output, dataset)
        for sample in os.listdir(speciesFolder):
            for gtee in gteeList:
                for rep in range(1, nReplicates+1):
                    speciesFile = os.path.join(speciesFolder,
                                                sample, gtee, 'R'+str(rep),
                                                args.speciesname)
                    geneFile = os.path.join(geneFolder, sample, gtee, 'R'+str(rep),
                                                args.genename)
                    outputFile = os.path.join(outputFolder, sample,
                                                gtee, 'R'+str(rep),
                                                args.outputname)
                    outputParent = Path(outputFile).parent
                    outputParent.mkdir(parents=True, exist_ok=True)

                    cmd = f'python {args.qpath} -t {speciesFile} -g {geneFile} -o {outputFile} 1>/dev/null'

                    retValue = os.system(cmd)
                    if retValue != 0:
                        print(f'Error running the following command:\n{cmd}')

