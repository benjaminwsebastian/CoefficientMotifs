#!/bin/python
import numpy as np
import argparse

###########
## INPUT ##
###########

def getArgs():
    parser = argparse.ArgumentParser(description="create data of degree n", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--n",type=int,required=True,help="create polynomial of highest degree n (x^1+...+x^n)")
    parser.add_argument("--l",type=int,required=False,help="number of data points",default=600)
    
    args = parser.parse_args()
    return args

##########
## MAIN ##
##########

args = getArgs()
n = args.n
print n
length = args.l

x = np.linspace(-5,5,length)
coeffs = np.random.uniform(-20,20,n)
y = np.poly1d(coeffs)(x)+np.random.uniform(-20,20,length)

OUT = open('n'+str(n)+'.txt','w')
for i in range(len(x)):
    OUT.write('%f\t%f\n' % (x[i],y[i]))
