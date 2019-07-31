#!/bin/python3
import string,json,argparse
import numpy as np
from bs4 import BeautifulSoup

'''
This script reads the meme output file and prints the motif as a sequence of values (within some error eps=max{len(Bin[0]),len(bin[len(bins)-1])}.
'''

###########
## INPUT ##
###########

def get_args():
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
def symlog_shift(arr, shift=0):
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

def find_motif_info(fileName):
    with open(fileName,'r') as f:
        contents = f.read()

        soup = BeautifulSoup(contents, 'lxml')
        text = soup.head.text
        front_trimmed = text.rsplit('=')[1][1:] # Have to trim because meme didn't format the html output nicely
        back_trimmed = front_trimmed.rsplit(';')[0]
        info = json.loads(back_trimmed)
        
        motif = info['motifs'][0]['id']
        nsites = info['motifs'][0]['nsites']
    
        return motif,nsites
    
##########
## MAIN ##
##########

args = get_args()
file_name = args.f
motif = args.motif
unlog = args.unlog
bins_file = args.bins

alphabet = list(string.ascii_letters)
IN = open(bins_file)
bins = json.load(IN)

if file_name != '':
    motif,nsites = find_motif_info(file_name)

approx = []
for element in motif:
    index = alphabet.index(element)
    approx.append((bins[index]+bins[index+1])/2.0)


print('motif: %s\nNsites: %s' % (motif,nsites))

if unlog:
    print(unlog10(approx))
else:
    print(approx)
