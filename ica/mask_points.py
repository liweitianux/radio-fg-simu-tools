#!/usr/bin/env python

import numpy
import scipy
import mdp
import sys
import math
import exceptions
import img_utils
import pyfits

#r=2.56
#r=7.68

r=5.12

def main(input_names,threshold,output_name):
    img=pyfits.open(input_names[0])[0].data
    
    mask=numpy.ones(img.shape,'int')


    for input_name in input_names:
        print input_name
        img=pyfits.open(input_name)[0].data
        for i in range(0,img.shape[0]):
            for j in range(0,img.shape[1]):
                if img_utils.is_peak(img,i,j) and img[i,j]>=threshold:
                    for k in img_utils.make_round_region(img,i,j,r):
                        mask[k]=0


    print 'masked ratio=',(1-float(scipy.sum(mask))/img.size)*100,'%'
    pyfits.HDUList(pyfits.PrimaryHDU(mask)).writeto(output_name,clobber=True)

    

if __name__=='__main__':
    if not len(sys.argv)>=4:
        print "Usage:", sys.argv[0], "<input image> [input image2] [input image3] ... <threshold> <output mask name>"
        sys.exit(-1)
    main(sys.argv[1:-2],float(sys.argv[-2]),sys.argv[-1])
                
    
