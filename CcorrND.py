import sys, re, os.path
import pandas as pd
import numpy as np
from numpy.linalg import norm
from scipy.interpolate import interp1d
from scipy.stats import zscore
from operator import itemgetter
import ntpath
import random

__author__ = 'Rob Moseley'
'''

Input Data Format:
    - 3 bioreps per 12 time point samples
    - rows = genes
    - columns = time points

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

Usage:
    - python CcorrND.py <input_data> <spline_number> <max_time_delay> <maxcount>

'''
sys.path.insert(0,'/Users/rkd/Desktop/Ccorr/Network-Deconvolution-python')
from ND import ND

def cxcorr(x, y):
    """
    For calculating the circular cross correlation between two vectors
    :param x: real or complex vector
    :param y: real or complex vector
    :return: circular correlation coefficients and respective lags
    """
    x = x/norm(x) # normalization
    y = y/norm(y) # normalization
    y = np.asarray(y)
    x = np.asarray(x)
    c = []
    for k in range(len(y)):
        cor = sum(x * y.conj().transpose())
        c.append(cor)
        # circular shift
        y = np.insert(y, 0, y[-1])
        y = y[:-1]
    x = range(len(y)) # lags
    return x, c

file = sys.argv[1] # <input_data>
# path = os.path.dirname(file)
# fileName = os.path.basename(file)[0]
fileName = os.path.splitext(ntpath.basename(file))[0]

if not os.path.exists(fileName + "output"):
    os.makedirs(fileName + "_output")

# get biorep data
data = pd.read_csv(file)
data = data.set_index("genes")

# copy last three reps to front
data = pd.concat([data.iloc[:, -3:], data], axis=1)

# average data
data = data.transpose()
data = data.groupby(np.arange(len(data))//3).mean()
data = data.transpose()
data.columns = ['Leaf_T06p1', 'Leaf_T08', 'Leaf_T10', 'Leaf_T12',
                'Leaf_T14', 'Leaf_T16', 'Leaf_T18', 'Leaf_T20',
                'Leaf_T22', 'Leaf_T00', 'Leaf_T02', 'Leaf_T04', 'Leaf_T06p2']
data.to_csv(fileName + "_output/" + fileName + "_AveragedData.txt", sep='\t')

# cubic spline interpolation
splineMatrix = []
numbTP = len(data.columns)
x = np.arange(numbTP)+1
splineNumb = int(sys.argv[2]) # <spline_number>
for idx, row in data.iterrows():
    xx = np.linspace(1, numbTP, splineNumb)
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
final_network = []
final_timeDelays = []
T = int(splineNumb) - 1
r = int(sys.argv[3]) # <max_time_delay>

for i in range(numbGenes):
    print str(i+1) + ". Determing potential regulators for: " + geneNames[i]
    Li = []
    for j in range(numbGenes):
        # print i, j
        lag, ccorr = cxcorr(zscoreMatrix.iloc[j,:],zscoreMatrix.iloc[i,:])
        # print "\t" + geneNames[j]   +": lag = " + str(lag) + " - ccorr = " + str(ccorr)
        # get max absolute ccorr value and respective index within max time delay
        maxAbsCcorr = max(ccorr[:r], key=abs)
        maxAbsCcorrIdx = np.argmax(np.absolute(ccorr[:r]))
        # get lag value
        Li.append(lag[maxAbsCcorrIdx])
    Xi = []
    # get vector of time samples for target gene (i)
    # print Li
    X = np.asarray(zscoreMatrix[i:i+1])
    X = np.concatenate((X[:,r:T], X[:,:r]),axis=1)
    Xi.append(X)
    # get vectors of time samples for potential regulators (j)
    for j in range(numbGenes):
        if j != i: # ignore target gene (i)
            Lij = Li[j]
            xj = np.asarray(zscoreMatrix[j:j+1])
            xj = np.concatenate((xj[:,r-Lij:], xj[:,:r-Lij]),axis=1) # extract time samples
            Xi.append(xj) # align time samples
    # stack sequences row wise and convert to dataframe
    Xi = pd.DataFrame(np.vstack(Xi)).transpose()
    # transpose so genes are columns and create spearman correlaton matrix (goes col by col)
    C = Xi.corr(method='spearman')
    # apply network deconvolution
    NetDe = ND(C.as_matrix())
    # get vector which contains the direct dependencies between i and its potential regulators
    # gene i will be first column/row
    NDi = np.asarray(NetDe[0,:])
    # adjust vector to align properly with original gene order
    NDi =  np.concatenate((NDi[numbGenes-i:], NDi[0:numbGenes-i]))
    # add to final network
    final_network.append(NDi)
    # add time delays to final matrix
    final_timeDelays.append(Li)

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
rankedList = {}
# extract ccorr and time delay for each gene pair
for i in range(len(FN)):
    for j in range(len(FN)):
        if final_network[i][j] != 0:
            rankedList[geneNames[i] + ',' + geneNames[j]] = (final_network[i][j],  final_timeDelays[i][j])
rankedList = sorted(rankedList.iteritems(), key=itemgetter(1), reverse=True)

maxCount = sys.argv[4] # <maxcount>
nWrite = len(rankedList)
if nWrite < maxCount > 0:
    nWrite = int(maxCount)

outfile = open(fileName + "_output/" + fileName + "_Ranked_List.txt", 'w')
outfile.write("Source\tTarget\tCcorr\tTimeDelay\n")
for i in range(nWrite):
    (TF, target) = re.split(r',', rankedList[i][0].rstrip())
    (CC, TD) = rankedList[i][1]
    CC = float(CC)
    TD = float(TD)
    outfile.write("%s\t%s\t%.5f\t%d\n" % (TF, target, CC, TD))
outfile.close()
