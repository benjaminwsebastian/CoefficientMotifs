#!/bin/python3
import numpy as np
import argparse

###########
## INPUT ##
###########

def get_args():
    parser = argparse.ArgumentParser(description="create data of degree n", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--n",type=int,required=True,help="create polynomial of highest degree n (x^1+...+x^n)")
    parser.add_argument("--l",type=int,required=False,help="number of data points",default=600)
    
    args = parser.parse_args()
    return args

##########
## MAIN ##
##########

args = get_args()
N = args.n
print(N)
LENGTH = args.l

x = np.linspace(-5,5,LENGTH)
coeffs = np.random.uniform(-20,20,N+1) # Define a random polynomial of highest degree n
y = np.poly1d(coeffs)(x)+np.random.uniform(-20,20,LENGTH) # Get polynomial points and add uniform noise

OUT = open('n'+str(N)+'.txt','w') # Print the data
for i in range(len(x)):
    OUT.write('%f\t%f\n' % (x[i],y[i]))
