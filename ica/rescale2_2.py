#!/usr/bin/env python

import numpy
import pyfits

for i in open('cluster_list.txt'):
    freq,clu_num=i.split()
  
    img1=pyfits.open('ica_'+freq+'_'+clu_num+'.fits')[0].data
    mask_excl=pyfits.open('mask_excl.fits')[0].data
#    mask_edge=pyfits.open('mask_edge.fits')[0].data
    
    import utils
    arr1=utils.two_to_one(img1,mask_excl)
    
    #img3=pyfits.open('model_'+freq+'.fits')[0].data
    img3=pyfits.open('/home/tsubasa/work/simu_diffuse/free/2.0sample/Timg_bin_'+freq+'_sm.fits')[0].data
    arr3=utils.two_to_one(img3,mask_excl)


    sx2=(arr1**2).sum()
    sxy=(arr1*arr3).sum()
    sx=(arr1).sum()
    sy=(arr3).sum()

    n=len(arr1)

    k=(sxy*n-sy*sx)/(sx2*n-sx**2)
    b=(sx2*sy-sx*sxy)/(sx2*n-sx**2)

    img1=k*img1+b
    print freq
    pyfits.HDUList(pyfits.PrimaryHDU(img1)).writeto(freq+'.fits',clobber=True) 

