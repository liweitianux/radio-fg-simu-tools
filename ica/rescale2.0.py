#!/usr/bin/env python

import numpy
import pyfits

for i in open('cluster_list.txt'):
    freq,clu_num=i.split()
  
    img1=pyfits.open('ica_'+freq+'_'+clu_num+'.fits')[0].data
    mask_excl=pyfits.open('mask_excl.fits')[0].data
    mask_edge=pyfits.open('mask_edge.fits')[0].data
    
    import utils
    arr1=utils.two_to_one(img1,mask_excl)
    print numpy.mean(arr1)
    
#img3=pyfits.open('../components_65/ptr_65_sm.fits')[0].data
    img3=pyfits.open('/home/tsubasa/work/simu_diffuse/free/2.0sample/Timg_bin_'+freq+'_sm.fits')[0].data
    #img3=pyfits.open('/home/tsubasa/work/simu_diffuse/cluster/sim2_1/Timg_cluABbin_'+freq+'_sm.fits')[0].data
    #img3=pyfits.open('model_'+freq+'.fits')[0].data
    arr3=utils.two_to_one(img3,mask_excl)

#    if freq=='95':

#        for j in range(0,len(arr3)):
#            print arr3[j],arr1[j]

    # XXX: Minus mean value by 2 times??
    arr1-=numpy.mean(arr1)
    arr1-=numpy.mean(arr1)
    arr3-=numpy.mean(arr3)
#img1=utils.one_to_two(arr1,mask_excl)
    
    
    
    print freq,(numpy.mean(arr3*arr1)/numpy.mean(arr1*arr1))
    
    img1*=(numpy.mean(arr3*arr1)/numpy.mean(arr1*arr1))
    img1-=numpy.max(utils.two_to_one(img1,mask_edge))

#    img1=(img1-img1.min())
#    print "max=",img1.max()

#    img1/=img1.max()
#    img1=img1*(max(arr3)-min(arr3))+min(arr3)
#    print "max=",img1.max()
    #img1-=numpy.max(utils.two_to_one(img1,mask_edge))

    print arr1.std(),arr3.std()

    pyfits.HDUList(pyfits.PrimaryHDU(img1)).writeto(freq+'.fits',clobber=True) 

