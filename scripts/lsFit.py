#!/bin/python
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import string,json,argparse

# TODO: put the poly in coeffs > poly[x] to coeffs[poly][x] OR change coeffs into a list -- see if doing that affects other parts of the code

'''
This script takes an tab delimited file of 2D data.

It does three things:

(1) Fits polynomials, from c*x to C*SUM(x^n) (1 <= n < infinity and C is a lower diagonal coefficient matrix) for a given max n, to the 2D data.

(2) Creates a histogram of the coefficients, a heatmap of the coefficients, and overlays plots of the fitted polynomials along with the original data. (The axis of the overlay plot are scaled using symlog as Runge's phenomenon is observed in overfitting leading to an exponentially increasing range.

(3) Bins the coefficients into 26 bins then assigns each bin a letter in A to Z, with each successive bin being the next letter in the alphabet. So the sequence of coefficients is turned into a sequence of letters. These are then printed in fasta format so the program "meme" (part of the MEME Suite, a motif-based sequence analysis tool which can detect denovo motifs- specifically made for DNA or RNA sequences [but a different alphabet can be supplied, which is what I do]) can be run on it.
'''
###########
## INPUT ##
###########

def getArgs():
    parser = argparse.ArgumentParser(description="translate motif back into numbers, can scrape meme.hmtl for motif or take a motif", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--f",type=str,required=True,help="input data file (tab delimited)",default='')
    parser.add_argument("--o",type=str,required=True,help="output name")
    parser.add_argument("--n",type=int,required=False,help="Optional: largest polynomial approx (x^1 to x^n)",default=100)
    parser.add_argument("--symlog",type=bool,required=False,help="Optional: symlog the data",default=False)
    parser.add_argument("--label",type=bool,required=False,help="Optional: label heatmap rows",default=True)
    
    args = parser.parse_args()
    return args

def readDataFile(dataFile):
    xs,ys = [],[]
    for line in open(dataFile,'r'):
        x,y = line.strip().split('\t')
        xs.append(float(x))
        ys.append(float(y))
    return xs,ys

#################
## SUBROUTINES ##
#################

'''
from matplotlib source code
'''
def symlogShift(arr, shift=0):
    # shift array-like to symlog array with shift
    logv = np.abs(arr)*(10.**shift)
    logv[np.where(logv<1.)] = 1.
    logv = np.sign(arr)*np.log10(logv)
    return logv 

##########
## MAIN ##
##########

args = getArgs()
dataFile = args.f
n0 = args.n
symlog = args.symlog
outName = args.o
label = args.label

highestPoly = n0

errors = []
polys = {}
coeffs = []

x,y = readDataFile(dataFile)

'''
Fit the polynomials to the data
'''
for i in range(highestPoly):
    coeffs.append(np.polyfit(x,y,i+1))
    polys[i] = np.poly1d(coeffs[i])

'''
Plot overlayed figure
'''
plt.figure()
plt.subplot(221)
for i in polys:
    plt.plot(x,polys[i](x),'b-',alpha=0.1)

plt.plot(x,y,'go',markersize = 1)

if symlog:
    plt.yscale('symlog')
plt.title('%s %d' % ('Poly Fit With Highest Degree',highestPoly))
plt.xlabel('x')
plt.ylabel('y')

'''
Plot heatmap of coefficients
'''
chart = []
for poly in coeffs:
    if symlog:
        coeffs[coeffs.index(poly)] = symlogShift(poly)
    chart.append( np.concatenate([ poly[::-1] , np.zeros(highestPoly-len(poly)+1) ]) ) # add coefficients in reverve order so they are smallest to largest with respect to their powers and zero fill

plt.subplot(212)
sns.heatmap(chart,yticklabels=False,annot_kws={"size":2}) # can add vmin= and vmax= for labeling the bar, center=0 for making the bar centered at 0, cmap='Blues' for an alternative color scheme, can define rowLabels and set yticklabels=rowLabels but not advised because it gets crowded quickly
plt.xlabel("coefficients")
plt.ylabel("highest degree (n)")
plt.title("Heat Map of the Size of Coefficients For Each Polynomial up to "+str(n0))

'''
Plot distribution of coefficients
''' 
coeffsList = np.array([])
for polyIndex in range(len(coeffs)):
    coeffsList = np.concatenate([ coeffsList, coeffs[polyIndex][0:polyIndex] ]) # only adding to the ith term because the remaining is zero-filled and keeping those would skew the histogram.
    
plt.subplot(222)
maxCoeff = max(coeffsList)
minCoeff = min(coeffsList)
binWidth = (maxCoeff-minCoeff)/26.

n, bins, bars = plt.hist(coeffsList,bins=np.arange(minCoeff,maxCoeff,binWidth),normed=1,color='r')
plt.title("Distribution of coefficients")
plt.xlabel("symlog of coefficient")
plt.tight_layout()
plt.savefig(outName+'_n'+str(n0)+'.pdf')

'''
Transform coeff data into a sequence of letters then print to fasta file
'''
alphabet=string.ascii_letters
seqs = []
for i in chart: # i is the sequence of coefficients corresponding to one fitted polynomial
    seq = []
    for j in range(len(i)):
        for k in range(len(bins)-1): # Loop through alphabet to determine the bin for coefficient "chart[i][j]"
            if i[j] >= bins[k] and i[j] <= bins[k+1]:
                seq.append(alphabet[k])
    seq = ''.join(seq)
    seqs.append(seq)

OUT = open(outName+'_n'+str(n0)+'_coeffSeqs.fasta','w')
for i in range(len(seqs)):
    OUT.write('>%d\n%s\n' % (i+1,seqs[i]))

OUT = open(outName+'_n'+str(n0)+'_bins.json','w')
json.dump(bins.tolist(),OUT)
