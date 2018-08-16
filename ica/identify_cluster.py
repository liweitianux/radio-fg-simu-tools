#!/usr/bin/env python

import math
import numpy
import pyfits
import sys

def kurt(m):
    m0=m-numpy.mean(m)
    E4=numpy.sum(m0**4)/m.size
    E2=numpy.sum(m0**2)/m.size

    return E4/E2/E2-3

if __name__=='__main__':
    if len(sys.argv)!=3:
        print "Usage:",sys.argv[0]," <prefix> <freq>"
        sys.exit(-1)


    max_kurt=0
    max_kurt_number=0
    prefix=sys.argv[1]
    freq=sys.argv[2]
    for i in range(1,5):
        fname=prefix+freq+'_'+`i`+'.fits'
        k=kurt(pyfits.open(fname)[0].data)

        if k>max_kurt:
            max_kurt=k
            max_kurt_number=i

    print '%s %s'%(freq,max_kurt_number)

    
        

