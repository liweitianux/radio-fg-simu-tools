#!/usr/bin/env python


import numpy
import scipy
import mdp
import sys
import math
import pyfits
import exceptions

class mask_icaer:
    source_list=[]
    mask=None
    
    def __init__(self,nout):
        self.nout=nout
        pass

    def two_to_one(self,mx):
        mask=self.mask
        size1=numpy.sum(mask)
        result=numpy.zeros([size1],mx.dtype)
        n=0
        for i in range(0,mask.shape[0]):
            for j in range(0,mask.shape[1]):
                if(mask[i,j]!=0):
                    result[n]=mx[i,j]
                    n+=1
        return result
    
                
    def one_to_two(self,mx):
        mask=self.mask
        size1=mx.shape[0]
        result=numpy.zeros(mask.shape,mx.dtype)
        n=0
#        print mx.shape
#        print size1
        for i in range(0,mask.shape[0]):
            for j in range(0,mask.shape[1]):
                if mask[i,j]!=0:
                    result[i,j]=mx[n]
                    n+=1
        return result


    
    def append_source(self,img):
        if self.mask.shape!=img.shape:
            raise exceptions.RuntimeError("shape do not matches")
        self.source_list.append(self.two_to_one(img))

    def perform(self):
        self.source_list=numpy.array(self.source_list,'double')
        self.source_list=self.source_list.transpose()
       
        guess_matrix=None
        self.source_list=mdp.fastica(self.source_list, approach='defl',g='gaus', guess=guess_matrix,fine_g='gaus', mu=0.01, stabilization=0.001, sample_size=1, fine_tanh=1, max_it=50000000, max_it_fine=100, failures=500, limit=0.001, verbose=False, whitened=False, white_comp=self.nout, white_parm={'svd':True})

        result=[]
        print self.source_list.shape
        for i in range(0,nout):
            result.append(self.one_to_two(self.source_list[:,i]))

        #sources=self.source_list
        self.source_list=[]
        return result


if __name__=="__main__":

    nsources=len(sys.argv)-4
    if nsources<2:
        print "Usage:",sys.argv[0]," <mask> <nout> <signal 1> <signal 2> [signal 3,4...] <out prefix>"

        sys.exit(-1)


    mask=pyfits.open(sys.argv[1])[0].data
    
    nout=int(sys.argv[2])
    out_prefix=sys.argv[3+nsources]
    print "out prefix:",out_prefix
    mi=mask_icaer(nout)
    mi.mask=(mask.astype('int'))
    for i in range(3,3+nsources):
        print "loading:",sys.argv[i]
        mi.append_source(pyfits.open(sys.argv[i])[0].data)

    print "Performing ICA"
    
    sources=mi.perform()

    for i in range(0,nout):
        outname=out_prefix+`i+1`+".fits"
        print "Writing ",outname
        pyfits.HDUList([pyfits.PrimaryHDU(sources[i])]).writeto(outname,clobber=True)


        
