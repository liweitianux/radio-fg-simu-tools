#!/usr/bin/env python

import numpy
import pyfits
import sys


def split_mean(arr):
    m=numpy.mean(arr)
    return arr-m,m
def usage():
    print "rescale -m <mask> [mask2] ... -s <source1> [source2] [...] -t <target1> [target2] [...] -o <output prefix>"



def two_to_one(mx,mask):
    size1=numpy.sum(mask)
    result=numpy.zeros([size1],mx.dtype)
    n=0
    for i in range(0,mask.shape[0]):
        for j in range(0,mask.shape[1]):
            if(mask[i,j]!=0):
                result[n]=mx[i,j]
                n+=1
    return result


def one_to_two(mx,mask):
    size1=mx.shape[0]
    result=numpy.zeros(mask.shape,mx.dtype)
    n=0
    for i in range(0,mask.shape[0]):
        for j in range(0,mask.shape[1]):
            if mask[i,j]!=0:
                result[i,j]=mx[n]
                n+=1
    return result


def calc_pp(mx,target):
    mx-=numpy.mean(mx)
    target-=numpy.mean(target)
    return numpy.mean(target*mx)/numpy.mean(mx*mx)



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
        image_size=rawfits.size
        image_shape=rawfits.shape
        arr=two_to_one(rawfits,mask1)
        target_mean.append(numpy.mean(arr))
        arr-=target_mean[-1]
        target_data.append(arr)
        
        
        
    for j in source_names:
        rawfits=pyfits.open(j)[0].data
        image_size=rawfits.size
        image_shape=rawfits.shape
        arr=two_to_one(rawfits,mask1)
        source_mean.append(numpy.mean(arr))
        arr-=source_mean[-1]
        source_data.append(arr)


    if len(source_mean)>len(target_mean):
        print "Warning: absolute level can not be determined!"
        exit()

    if len(source_mean)==0:
        usage()
        exit()

    if outprefix==None:
        usage()
        exit()


    pp_matrix=numpy.zeros([len(target_mean),len(source_mean)])
    
    for i in range(0,len(target_data)):
        for j in range(0,len(source_data)):
            pp_matrix[i,j]=calc_pp(source_data[j],target_data[i])

    print pp_matrix
    [l,resd,rank,s]=numpy.linalg.lstsq(pp_matrix,target_mean)
    print l
    for i in range(0,len(source_names)):
        outname="%s%d.fits"%(outprefix,i+1)
#        print source_names[i],' ',pp_matrix[

        img=pyfits.open(source_names[i])[0].data
        img+=l[i]
        img*=pp_matrix[0,i]
        
        #img*=mask1
        pyfits.HDUList([pyfits.PrimaryHDU(img)]).writeto(outname,clobber=True)
    
if __name__=='__main__':
    main(sys.argv)
