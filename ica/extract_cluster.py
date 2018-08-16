#!/usr/bin/env python

import re
import numpy
import scipy
import mdp
import sys
import math
import exceptions
import img_utils
import pyfits


#mask=numpy.zeros(pyfits.open(sys.argv[2])[0].data.shape,'int')

if len(sys.argv)!=4:
    print "Usage:",sys.argv[0], "<region> <input> <output>"
    sys.exit(-1)


line=open(sys.argv[1]).readlines()[0]
splited=re.findall(r"(\S+)\((\S+),(\S+),(\S+)\)",line.strip())
print splited
center_j=float(splited[0][1])
center_i=float(splited[0][2])
r=float(splited[0][3])

img=pyfits.open(sys.argv[2])[0].data

outimg=numpy.zeros([r*2+1,r*2+1])
outimg*=float('nan')
for di in range(-r,r+1):
    for dj in range(-r,r+1):
        i=di+center_i
        j=dj+center_j

        if di*di+dj*dj<r*r:
            outimg[di+r,dj+r]=img[i,j]
        else:
            outimg[di+r,dj+r]=float("nan")

pyfits.HDUList(pyfits.PrimaryHDU(outimg)).writeto(sys.argv[3],clobber=True)
