import matplotlib
matplotlib.use('pgf')
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{amsfonts}\usepackage{amssymb}')
import numpy as np
import pickle

'''
Change the filenames to the pickle files saved by computeCDStats.py
'''
if __name__ == '__main__':
    bestDataFileName = '../junk/best-rooting-astral.pickle'
    with open(bestDataFileName, 'rb') as dataFile:
        bestOutputData = pickle.load(dataFile)

    qrleDataFileName = '../junk/qr-le.pickle'
    with open(qrleDataFileName, 'rb') as qrleFile:
        qrleData = pickle.load(qrleFile)

    qruDataFileName = '../junk/qru.pickle'
    with open(qruDataFileName, 'rb') as qruFile:
        qruData = pickle.load(qruFile)

    datasets = [str(i)+'_taxon' for i in [5, 10, 15]]
    gteeList = ['true', '1500', '1000', '500', '250']
    xList = np.array([0.5, 1.5, 2.5, 3.5, 4.5])
    yhighs = [0.65, 0.65, 0.65]
    for datasetNum,dataset in enumerate(datasets):
        plt.figure()
        try:
            meanVal = np.zeros(len(gteeList))
            errorVal = np.zeros(len(gteeList))
            for i,gtee in enumerate(gteeList):
                meanVal[i] = np.mean(bestOutputData[dataset][gtee])
                errorVal[i] = np.std(bestOutputData[dataset][gtee])
            plt.errorbar(xList-0.2, meanVal, yerr=errorVal, fmt='o', color='black',
                            ecolor='lightgray', elinewidth=3, capsize=0, label='best rooting of ASTRAL')
        except:
            pass
        try:
            meanVal = np.zeros(len(gteeList))
            errorVal = np.zeros(len(gteeList))
            for i,gtee in enumerate(gteeList):
                meanVal[i] = np.mean(qrleData[dataset][gtee])
                errorVal[i] = np.std(qrleData[dataset][gtee])
            plt.errorbar(xList+0.2, meanVal, yerr=errorVal, fmt='o', color='red',
                            ecolor='lightgray', elinewidth=3, capsize=0, label='QR (with ASTRAL trees)')
        except:
            pass
        try:
            meanVal = np.zeros(len(gteeList))
            errorVal = np.zeros(len(gteeList))
            for i,gtee in enumerate(gteeList):
                meanVal[i] = np.mean(qruData[dataset][gtee])
                errorVal[i] = np.std(qruData[dataset][gtee])
            plt.errorbar(xList, meanVal, yerr=errorVal, fmt='o', color='blue',
                            ecolor='lightgray', elinewidth=3, capsize=0, label='QRU (with ASTRAL trees)')
        except:
            pass
        #plt.grid()
        plt.title(f'{dataset[:-6]}-taxon data', fontsize='x-large')
        tickLabels = ['true']
        for gtee in gteeList[1:]:
            tickLabels.append(f'$L={gtee}$')
        dividers = 0.5*(xList[1:] + xList[:-1])
        for val in dividers:
            plt.plot([val, val], [-0.1, np.max(yhighs)+0.1], 'k--', linewidth=0.5)
        plt.xticks(ticks=xList, labels=tickLabels)
        plt.xlim([0., 5.])
        plt.ylim([0., yhighs[datasetNum]])
        plt.tick_params(top=False, bottom=False, left=True, right=False, labelsize='x-large')
        plt.ylabel(r'Clade distance', fontsize='large')
        plt.xlabel(r'(Length of MSA for gene-tree estimation)', fontsize='large')
        plt.legend(fontsize='large', loc='upper left', framealpha=1.0)
        plt.savefig(f'{dataset}-plot.pdf')

