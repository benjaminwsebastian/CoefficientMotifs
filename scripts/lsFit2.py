#!/bin/python
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import string,json,argparse

###########
## INPUT ##
###########

def getArgs():
    parser = argparse.ArgumentParser(description="translate motif back into numbers, can scrape meme.hmtl for motif or take a motif", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--f",type=str,required=False,help="input data file (tab delimited)",default='')
    parser.add_argument("--n",type=int,required=True,help="largest polynomial approx (x^1 to x^n)")
    parser.add_argument("--symlog",type=bool,required=False,help="symlog the data",default=False)
    parser.add_argument("--o",type=str,required=True,help="output name")
    parser.add_argument("--label",type=bool,required=False,help="label heatmap rows",default=True)
    
    args = parser.parse_args()
    return args

def readDataFile(dataFile):
    xs,ys = [],[]
    for line in open(dataFile,'r'):
        x,y,z = line.strip().split('\t')
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
coeffs = {}

x,y = readDataFile(dataFile)

for i in range(highestPoly):
    coeffs[i] = np.polyfit(x,y,i+1)
    polys[i] = np.poly1d(coeffs[i])
    errors.append( (i,np.sqrt(np.mean( (y-polys[i](x))**2 ))) )

errors = sorted(errors,key=lambda p:p[1])

bestFit = errors[0][0]

'''
Plot overlayed figure
'''
plt.figure()
plt.subplot(221)
for i in polys:
    plt.plot(x,polys[i](x),'b-',alpha=0.1)

plt.plot(x,y,'go',markersize = 1)
plt.plot(x,polys[bestFit](x),'r-')

if symlog:
    plt.yscale('symlog')
plt.title('%s %d \n%s %d %s %f' % ('Best Poly Fit With Highest Degree',highestPoly,'Best fit at n=',bestFit+1,'and error ',errors[0][1]))
plt.xlabel('x')
plt.ylabel('y')

'''
Plot heatmap of coefficients
'''
rowLabels = []
if label:
    for poly in errors:
        rowLabels.append(str(poly[0]+1))
else:
    rowLabels = False
        
chart = []
for poly in errors:
    if symlog:
        coeffs[poly[0]] = symlogShift(coeffs[poly[0]])
    chart.append( np.concatenate([ coeffs[poly[0]] , np.zeros(highestPoly-poly[0]-1) ]) ) # add coefficients and zero fill

plt.subplot(212)
sns.heatmap(chart,yticklabels=rowLabels,) # can add vmin= and vmax= for labeling the bar and center=0 cmap='Blues'
plt.xlabel("coefficients")
plt.ylabel("highest degree (n)")
plt.title("Heat Map of the Size of Coefficients For Each Polynomial up to "+str(n0))

'''
Plot distribution of coefficients
''' 
coeffsList = np.array([])
for i in coeffs:
    coeffsList = np.concatenate([ coeffsList, coeffs[i] ])
    
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
Transform coeff data into letters then print to fasta file so I can run meme on it
'''

alphabet=string.ascii_letters
seqs = []
for i in chart:
    seq = []
    for j in range(len(i)):
        for k in range(len(bins)-1):
            if i[j] >= bins[k] and i[j] <= bins[k+1]:
                seq.append(alphabet[k])
    seq = ''.join(seq)
    seqs.append(seq)

OUT = open(outName+'_n'+str(n0)+'_coeffSeqs.fasta','w')
for i in range(len(seqs)):
    OUT.write('>%d\n%s\n' % (errors[i][0]+1,seqs[i]))

OUT = open(outName+'_n'+str(n0)+'_bins.json','w')
json.dump(bins.tolist(),OUT)
