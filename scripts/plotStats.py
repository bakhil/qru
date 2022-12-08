import matplotlib
matplotlib.use('pgf')
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{amsfonts}\usepackage{amssymb}')
import numpy as np
import scipy.stats as sps
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
    legendPltsList, legendStrList = [], []
    boxWidth = 0.16
    for datasetNum,dataset in enumerate(datasets):
        fig, ax = plt.subplots()
        legendPltsList, legendStrList = [], []
        def add_plot(data, shift, boxcolor, meancolor, legendStr):
            try:
                retPlt = ax.boxplot([data[gtee] for gtee in gteeList],
                                    positions=xList+shift, showmeans=True, meanline=True,
                                    patch_artist=True, widths=boxWidth,
                                    boxprops=dict(facecolor=boxcolor, color=meancolor),
                                    meanprops=dict(color=meancolor, linestyle='-'),
                                    medianprops=dict(linewidth=0))
                legendPltsList.append(retPlt['boxes'][0])
                legendStrList.append(legendStr)
            except Exception as e:
                print(e)
            return
        add_plot(bestOutputData[dataset], -0.2, 'lightgray', 'black', 'best rooting of ASTRAL')
        add_plot(qrleData[dataset], 0.2, 'cyan', 'navy', 'QR (with ASTRAL trees)')
        add_plot(qruData[dataset], 0., 'tan', 'darkred', 'QRU (with ASTRAL trees)')

        ax.set_title(f'{dataset[:-6]}-taxon data', fontsize='x-large')
        tickLabels = ['true']
        for gtee in gteeList[1:]:
            tickLabels.append(f'$L={gtee}$')
        dividers = 0.5*(xList[1:] + xList[:-1])
        for val in dividers:
            ax.plot([val, val], [-0.1, np.max(yhighs)+0.1], 'k--', linewidth=0.5)
        ax.legend(legendPltsList, legendStrList, fontsize='large', loc='upper left', framealpha=1.0)
        ax.set_xticks(ticks=xList)
        ax.set_xticklabels(labels=tickLabels)
        ax.set_xlim([0., 5.])
        ax.set_ylim([0., yhighs[datasetNum]])
        ax.tick_params(top=False, bottom=False, left=True, right=False, labelsize='x-large')
        ax.set_ylabel(r'Clade distance', fontsize='large')
        ax.set_xlabel(r'(Length of MSA for gene-tree estimation)', fontsize='large')
        plt.savefig(f'{dataset}-plot.pdf')

