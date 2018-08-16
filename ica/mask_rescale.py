#!/usr/bin/env python


import getopt
from scipy import *
import pyfits
import numpy
import os
import mdp
import sys


def split_mean(arr):
    m=numpy.mean(arr)
    return arr-m,m
def usage():
    print "rescale -m <mask> [mask2] ... -s <source1> [source2] [...] -t <target1> [target2] [...] -o <output prefix>"





def two_to_one(mx,mask):
    size1=sum(mask)
    result=zeros([size1],mx.dtype)
    n=0
    for i in range(0,mask.shape[0]):
        for j in range(0,mask.shape[1]):
            if(mask[i,j]!=0):
                result[n]=mx[i,j]
                n+=1
    return result


def one_to_two(mx,mask):
    size1=mx.shape[0]
    result=zeros(mask.shape,mx.dtype)
    n=0
    for i in range(0,mask.shape[0]):
        for j in range(0,mask.shape[1]):
            if mask[i,j]!=0:
                result[i,j]=mx[n]
                n+=1
    return result




def main(args):
    image_shape=None
    image_size=0
    source_data=[]
    target_data=[]
    source_mean=[]
    target_mean=[]

    source_names=[]
    target_names=[]
    mask_names=[]
    

    outprefix=[]

    current_ptr=None

    for i in args:
        if i=='-m':
            current_ptr=mask_names
            continue
        elif i=='-s':
            current_ptr=source_names
            continue
        elif i=='-t':
            current_ptr=target_names
            continue
        elif i=='-o':
            current_ptr=outprefix
            continue

        if current_ptr!=None:
            current_ptr.append(i)

    print "sources:"
    for i in source_names:
        print i
        
    print "targets:"
    for i in target_names:
        print i
        
    print "masks:"
    for i in mask_names:
        print i

    if len(source_names)<1:
        usage()
        sys.exit(-1)

    if len(target_names)<1:
        usage()
        sys.exit(-1)

    if len(mask_names)<1:
        usage()
        sys.exit(-1)

    if len(outprefix)!=1:
        usage()
        sys.exit(-1)

    outprefix=outprefix[0]
    
    if image_shape is None:
        image_shape=pyfits.open(mask_names[0])[0].data.shape

    mask1=numpy.ones(image_shape,'int')
    print mask1.shape
    for i in mask_names:
        mask2=pyfits.open(i)[0].data.astype('int')
        print mask2.shape
        mask1*=mask2
    


    for j in target_names:
        rawfits=pyfits.open(j)[0].data
        [rawfits,m]=split_mean(rawfits)
        target_mean.append(m)
        image_size=rawfits.size
        image_shape=rawfits.shape
        target_data.append(two_to_one(rawfits,mask1))
                
        
    for j in source_names:
        rawfits=pyfits.open(j)[0].data
        [rawfits,m]=split_mean(rawfits)
        source_mean.append(m)
        image_size=rawfits.size
        image_shape=rawfits.shape
        source_data.append(two_to_one(rawfits,mask1))


    if len(source_mean)>len(target_mean):
        print "Warning: absolute level can not be determined!"
        exit()

    if len(source_mean)==0:
        usage()
        exit()

    if outprefix==None:
        usage()
        exit()

    target_data=array(target_data).transpose()
    source_data=array(source_data).transpose()

    [w,resd,rank,s]=numpy.linalg.lstsq(source_data,target_data)

    print w
#    print w.shape
    print "target mean:",target_mean
    [l,resd,rank,s]=numpy.linalg.lstsq(w.transpose(),target_mean)
    print "source mean:",l
    print "factors:",w[:,0]
#    print sum(l*w[:,1])
    for i in range(1,len(source_mean)+1):
        outname="%s%d.fits"%(outprefix,i)
        print outname
#        img=one_to_two(source_data[:,i-1]+l[i-1])
        img=pyfits.open(source_names[i-1])[0].data
        img+=l[i-1]
        
        img*=w[i-1,0]
        img*=mask1
        pyfits.HDUList([pyfits.PrimaryHDU(img)]).writeto(outname,clobber=True)
        
    
    pass



if __name__=="__main__":
    main(sys.argv)


    
