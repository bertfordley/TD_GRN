# def CcorrND(input_file, s, r, maxCount):
# uncomment above and indent all below to make function
# comment sys.argv lines below, also
'''

Usage:
    - python CcorrND.py input_data s r maxCount

Input Arguments:

:input_file         Data Format:3 bioreps per 12 time point samples - rows = genes; columns = time points
:s                  spline number
:r                  max time delay
:maxCount           maximum number of ranked edges to save to file

Steps:
    - computes averages for each time point sample
    - converts data to 1-hour intervals
    - converts data to z-scores
    - runs cross correlation and network deconvolution for each gene
    - ranks gene-gene interactions

Outputs:
    - Averaged data
    - Z-Scores
    - Direct dependencies matrix
    - Time Delay matrix
    - Ranked List of gene-gene interactions

'''

import sys, os.path
import pandas as pd
import numpy as np
from numpy.linalg import norm
from scipy.interpolate import interp1d
from scipy.stats import zscore
from operator import itemgetter
import ntpath

__author__ = 'Rob Moseley'

sys.path.insert(0, '/Users/rkd/Desktop/Ccorr/Network-Deconvolution-python')
from ND import ND


# sys.path.insert(0,'/Users/rkd/Desktop/Ccorr/Network-Deconvolution-python')
# from test_ND import ND

# comment lines out for making function
input_file = sys.argv[1]
s = int(sys.argv[2])
r = int(sys.argv[3])
maxCount = int(sys.argv[4])


def cxcorr(a, b):
    """
    For calculating the circular cross correlation between two vectors
    :param a: real or complex vector
    :param b: real or complex vector
    :return: circular correlation coefficients and respective lags
    """
    a /= norm(a)  # normalization
    b /= norm(b)  # normalization
    b = np.asarray(b)
    a = np.asarray(a)
    c = []
    for k in range(len(b)):
        cor = sum(a * b.conj().transpose())
        c.append(cor)
        # circular shift
        b = np.insert(b, 0, b[-1])
        b = b[:-1]
    a = range(len(b))  # lags
    return a, c


# path = os.path.dirname(file)
# fileName = os.path.basename(file)[0]
fileName = os.path.splitext(ntpath.basename(input_file))[0]

if not os.path.exists(fileName + "output"):
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

# print zscoreMatrix.corr(method='spearman')

# run algorithm on z-scores
numbGenes = len(zscoreMatrix.index)
geneNames = zscoreMatrix.index
final_network = []
final_timeDelays = []
T = int(s) - 1
idxMemory = range(zscoreMatrix.shape[0])

for i in range(numbGenes):
    print str(i + 1) + ". Determing potential regulators for: " + geneNames[i]
    Li = []
    maxCC = []
    for j in range(numbGenes):
        # print i, j
        lag, ccorr = cxcorr(zscoreMatrix.iloc[j, :], zscoreMatrix.iloc[i, :])
        # print "\t" + geneNames[j]   +": lag = " + str(lag) + " - ccorr = " + str(ccorr)
        # get max absolute ccorr value and respective index within max time delay
        maxAbsCcorr = max(ccorr[:r], key=abs)
        maxAbsCcorrIdx = np.argmax(np.absolute(ccorr[:r]))
        maxCC.append(abs(maxAbsCcorr))
        # get lag value
        Li.append(lag[maxAbsCcorrIdx])

    geneIdx = 0
    for k in range(numbGenes):
        if k == i:
            geneIdx = k
	# initialize matrix with zeros
    Xi = np.zeros(shape=(len(idxMemory), zscoreMatrix.shape[1]))
    # get vector of time samples for target gene (i)
    X = np.asarray(zscoreMatrix[geneIdx:geneIdx + 1])
    X = np.concatenate((X[:, r:T], X[:, :r]), axis=1)
    Xi[[idxMemory.index(x) for x in idxMemory if x == geneIdx], :] = X
    # get vectors of time samples for potential regulators (j)
    for j in range(len(idxMemory)):
        jIdx = idxMemory[j]
        if jIdx != geneIdx:  # ignore target gene (i)
            Lij = Li[jIdx]
            xj = np.asarray(zscoreMatrix[jIdx:jIdx + 1])
            xj = np.concatenate((xj[:, r - Lij:], xj[:, :r - Lij]), axis=1)  # extract time samples
            Xi[[idxMemory.index(x) for x in idxMemory if x == jIdx], :] = xj  # align time samples
    # stack sequences row wise and convert to dataframe
    Xi = pd.DataFrame(np.vstack(Xi)).transpose()
    print "Number of genes: " + str(Xi.shape[1])
    # add vectors of zeros if Xi only contains gene i
    # transpose so genes are columns and create spearman correlaton matrix (goes col by col)
    C = Xi.corr(method='spearman')
    C.to_csv(fileName + "_output/" + fileName + "_spearmanMat.txt", sep="\t")
    # apply network deconvolution
    NetDe = ND(C.as_matrix())
    # get vector which contains the direct dependencies between i and its potential regulators
    # gene i will be first column/row
    NDi = np.asarray(NetDe[i, :])
    # add to final network
    final_network.append(NDi)
    # add time delays to final matrix
    final_timeDelays.append(np.asarray(Li))

# Save network and time delays
FN = pd.DataFrame(final_network)
FN.index = geneNames
FN.columns = geneNames
FN.to_csv(fileName + "_output/" + fileName + "_GRN.txt", sep='\t')

FD = pd.DataFrame(final_timeDelays)
FD.index = geneNames
FD.columns = geneNames
FD.to_csv(fileName + "_output/" + fileName + "_TimeDelays.txt", sep='\t')

print "Ranking interactions..."
rankedList = []
# extract ccorr and time delay for each gene pair
for i in range(len(final_network)):
    for j in range(len(final_network)):
        if final_network[i][j] >= 0.7:
            rankedList.append([geneNames[i], geneNames[j], final_network[i][j], final_timeDelays[i][j]])
rankedList = sorted(rankedList, key=itemgetter(2), reverse=True)

nWrite = len(rankedList)
if nWrite < maxCount > 0:
    nWrite = int(maxCount)

outfile = open(fileName + "_output/" + fileName + "_Ranked_List.txt", 'w')
outfile.write("Source\tTarget\tCcorr\tTimeDelay\n")
for i in range(nWrite):
    (TF, target, CC, TD) = rankedList[i]
    CC = float(CC)
    TD = float(TD)
    outfile.write("%s\t%s\t%.5f\t%d\n" % (TF, target, CC, TD))
outfile.close()
