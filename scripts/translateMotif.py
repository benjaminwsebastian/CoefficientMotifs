#!/bin/python
import string,json,argparse
import numpy as np
from bs4 import BeautifulSoup

###########
## INPUT ##
###########

def getArgs():
    parser = argparse.ArgumentParser(description="translate motif back into numbers, can scrape meme.hmtl for motif or take a motif", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--f",type=str,required=False,help="meme output file (html)",default='')
    parser.add_argument("--motif",type=str,required=False,help="motif",default='')
    parser.add_argument("--unlog",type=bool,required=False,help="unlog data (it is symlogged)",default=False)
    parser.add_argument("--bins",type=str,required=True,help="bin json file")
    
    args = parser.parse_args()                                                 
    return args    

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

def unlog10(arr):
    unlogged = []
    for coeff in arr:
        sgn = np.sign(coeff)
        unlogged.append(sgn*(10**coeff))
    return unlogged

def findMotif(fileName):
    with open(fileName,'r') as f:
        contents = f.read()

        soup = BeautifulSoup(contents, 'lxml')
        text = soup.head.text
        frontTrimmed = text.rsplit('=')[1][1:] # Have to trim because meme didn't format the html output nicely
        backTrimmed = frontTrimmed.rsplit(';')[0]
        info = json.loads(backTrimmed)
        motif = info['motifs'][0]['id']

        return motif
    
##########
## MAIN ##
##########

args = getArgs()
fileName = args.f
motif = args.motif
unlog = args.unlog
binsFile = args.bins

alphabet = list(string.ascii_letters)
IN = open(binsFile)
bins = json.load(IN)

if fileName != '':
    motif = findMotif(fileName)

approx = []
for element in motif:
    index = alphabet.index(element)
    approx.append((bins[index]+bins[index+1])/2.0)


print 'motif:',motif

if unlog:
    print unlog10(approx)
else:
    print approx
