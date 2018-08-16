#!/usr/bin/env python

#used to compare the maps in different frequencys

import numpy
import scipy
import pyfits
import sys
import math

img1=pyfits.open(sys.argv[1])[0].data
img2=pyfits.open(sys.argv[2])[0].data
#print numpy.corrcoef(img1.resize(img1.size),img2.resize(img2.size))
sx2=(img1**2).sum()
sxy=(img1*img2).sum()
sx=(img1).sum()
sy=(img2).sum()

n=img1.size

k=(sxy*n-sy*sx)/(sx2*n-sx**2)
b=(sx2*sy-sx*sxy)/(sx2*n-sx**2)
#mean=max(img1.mean(),img2.mean())
mean=img2.mean()
img1=k*img1+b

#print k,b
rms=math.sqrt(((img1-img2)**2).mean())

rrms=math.sqrt(((img2-img1)**2).sum()/(img2**2).sum())
#print rrms

#print rms/mean
#print scipy.corrcoef(img1.resize(img1.size),img2.resize(img2.size))
img2.resize(img2.size)
img1.resize(img1.size)
cc=scipy.corrcoef(img1,img2)
print cc[0,1]
