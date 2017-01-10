
'''

Usage:
    - python CcorrND.py <input_file> <s> <r> <thres>

Input Arguments:

:input_file         Data Format:3 bioreps per 12 time point samples - rows = genes; columns = time points
:s                  max value for linearly spaced vector (for cubic spline interoplation)
:r                  max time delay
:thres              threshold for skipping genes in time sample alignment (for cxcorr values)

Steps:
    - computes averages for each time point sample
    - converts data to 1-hour intervals
    - converts data to z-scores
    - runs cross correlation to find time delays
    - align time samples and perform spearman rank then network deconvolution for each gene

Outputs:
    - Averaged data
    - Z-Scores
    - Direct dependencies matrix in json format

'''

import sys, os.path
import pandas as pd
import numpy as np
from numpy.linalg import norm
from numpy.fft import fft, ifft
from scipy.interpolate import interp1d
from scipy.stats import zscore
from itertools import izip
import ntpath
from scipy.stats import spearmanr
import json

__author__ = 'Rob Moseley'

input_file = sys.argv[1]
s = int(sys.argv[2])
r = int(sys.argv[3])
thres = float(sys.argv[4])

def cxcorr(a, b):
    """
    For calculating the circular cross correlation between two vectors
    :param a: real or complex vector
    :param b: real or complex vector
    :return: circular correlation coefficients and respective lags
    """
    a /= norm(a)  # normalization
    b /= norm(b)  # normalization
    coeff = ifft(fft(a) * fft(b).conj()).real
    lag = range(len(b))
    return lag, coeff

fileName = os.path.splitext(ntpath.basename(input_file))[0]

if not os.path.exists(fileName + "_output"):
    os.makedirs(fileName + "_output")

# get biorep data
data = pd.read_csv(input_file)
data = data.set_index("genes")

# copy last three reps to front
data = pd.concat([data.iloc[:, -3:], data], axis=1)

# average data
data = data.transpose()
data = data.groupby(np.arange(len(data)) // 3).mean()
data = data.transpose()
data.columns = ['Leaf_T06p1', 'Leaf_T08', 'Leaf_T10', 'Leaf_T12',
                'Leaf_T14', 'Leaf_T16', 'Leaf_T18', 'Leaf_T20',
                'Leaf_T22', 'Leaf_T00', 'Leaf_T02', 'Leaf_T04', 'Leaf_T06p2']
data.to_csv(fileName + "_output/" + fileName + "_AveragedData.txt", sep='\t')

# cubic spline interpolation
splineMatrix = []
numbTP = len(data.columns)
x = np.arange(numbTP) + 1
for idx, row in data.iterrows():
    xx = np.linspace(1, numbTP, s)
    cs = interp1d(x, row, kind='cubic')(xx)
    splineMatrix.append(cs)

splineMatrix = pd.DataFrame(splineMatrix)
# remove repeated time sample (6am) and leave new 7am time sample
splineMatrix.drop(splineMatrix.columns.values[0], axis=1, inplace=True)
splineMatrix.columns = ['7am', '8am', '9am', '10am', '11am', '12am',
                        '1pm', '2pm', '3pm', '4pm', '5pm', '6pm',
                        '7pm', '8pm', '9pm', '10pm', '11pm', '12pm',
                        '1am', '2am', '3am', '4am', '5am', '6am']
splineMatrix.index = data.index

# convert to z-scores
zscoreMatrix = splineMatrix.apply(zscore, axis=1)
zscoreMatrix.to_csv(fileName + "_output/" + fileName + "_zScores.txt", sep='\t')

# run algorithm on z-scores
numbGenes = len(zscoreMatrix.index)
geneNames = zscoreMatrix.index
T = int(s) - 1
corrdDict = {}

for i in range(numbGenes):
    colDict = {}
    print str(i + 1) + ". Determining potential regulators for: " + geneNames[i]
    Li = []
    maxCC = []
    for j in range(numbGenes):
        lag, ccorr = cxcorr(zscoreMatrix.iloc[j, :], zscoreMatrix.iloc[i, :])
        # get max absolute ccorr value and respective index within max time delay
        maxAbsCcorr = max(ccorr[:r + 1], key=abs)
        maxAbsCcorrIdx = np.argmax(np.absolute(ccorr[:r + 1]))
        maxCC.append(abs(maxAbsCcorr))
        # get lag value
        Li.append(lag[maxAbsCcorrIdx])
    print "Thresholding..."
    idxMemory = []
    adjLi = []
    geneIdx = 0
    for k in range(numbGenes):
        if maxCC[k] >= thres:
            if k == i:
                geneIdx = k
            idxMemory.append(k)
            adjLi.append(Li[k])
    print "Aligning time samples..."
    Xi = np.zeros(shape=(len(idxMemory), zscoreMatrix.shape[1]))
    # get vector of time samples for target gene (i)
    X = np.asarray(zscoreMatrix[geneIdx:geneIdx + 1])
    X = np.concatenate((X[:, r + 1:T], X[:, :r + 1]), axis=1)
    Xi[[idxMemory.index(x) for x in idxMemory if x == geneIdx], :] = X
    # get vectors of time samples for potential regulators (j)
    for j in range(len(idxMemory)):
        jIdx = idxMemory[j]
        if jIdx != geneIdx:  # ignore target gene (i)
            Lij = adjLi[j]
            xj = np.asarray(zscoreMatrix[jIdx:jIdx + 1])
            xj = np.concatenate((xj[:, r + 1 - Lij:], xj[:, :r + 1 - Lij]), axis=1)  # extract time samples
            Xi[[idxMemory.index(x) for x in idxMemory if x == jIdx], :] = xj  # align time samples
    # stack sequences row wise and convert to dataframe
    Xi = pd.DataFrame(np.vstack(Xi)).transpose()
    print "Number of genes: " + str(Xi.shape[1])
    # add vectors of zeros if Xi only contains gene i
    if Xi.shape[1] == 1:
        colDict[i] = (1, 0, 0)

    elif Xi.shape[1] == 2:
        # mess up here, self regulation numbers are wrong, but they will be dropped later on anyway
        C, pval = spearmanr(Xi[0], Xi[1])
        for col, td in izip(idxMemory, adjLi):
            colDict[col] = (C, pval, td)

    else:
        C, pval = spearmanr(Xi)
        ###########  apply ND to spearman matrix ###########
        Ci = C[idxMemory.index(i)]
        pvali = pval[idxMemory.index(i)]
        for col, s, p, td in izip(idxMemory, Ci, pvali, adjLi):
            colDict[col] = (s, p, td)

    corrdDict[i] = colDict

with open(fileName + "_output/" + fileName + "_output.json", "w") as outfile:
    json.dump(corrdDict, outfile, sort_keys=True, indent=2)
