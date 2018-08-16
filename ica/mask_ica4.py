#!/usr/bin/env python


import numpy
import scipy
import mdp
import sys
import math
import pyfits
import exceptions


indices=[]

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
        global indices
        guess_matrix=None
        try:
            matrix_file=open("ica_matrix.dat")
            guess_matrix=numpy.fromfile(matrix_file,sep=" ")
            guess_matrix=guess_matrix.reshape([self.nout,self.nout])
            print guess_matrix
#            print "guess matrix loaded"
#            print guess_matrix
            matrix_file.close()
        except IOError:
            pass

    
        
        self.source_list=numpy.array(self.source_list,'double')
        self.source_list=self.source_list.transpose()
        input_dim=self.source_list.shape[1]
        print "input dim:",self.source_list.shape
#        try:
#            matrix_file=open("whiten_matrix.dat")
#            whiten_matrix=numpy.fromfile(matrix_file,sep=" ")
#            whiten_matrix=whiten_matrix.reshape([input_dim,self.nout])
#            for i in range(0,input_dim):
#                self.source_list[:,i]-=numpy.mean(self.source_list[:,i])
#            whitened_source=mdp.utils.mult(self.source_list,whiten_matrix)
#            matrix_file.close()
#        except IOError:
        wd=mdp.nodes.WhiteningNode(output_dim=self.nout,svd=True)
        wd.train(self.source_list)
        whitened_source=wd.execute(self.source_list)
        whiten_matrix=wd.get_projmatrix()
#        whiten_matrix.tofile("whiten_matrix.dat",sep=" ")
#            print whitened_source.shape
            
        
        
        #self.source_list=mdp.fastica(whitened_source, approach='symm',g='tanh', guess=guess_matrix,fine_g='tanh', mu=0.0001, stabilization=0.0000001, sample_size=1, fine_tanh=1, max_it=1000000, max_it_fine=100000, failures=5000, limit=0.0000001, verbose=False, whitened=True)

        ficad=mdp.nodes.FastICANode(approach='symm',g='tanh', guess=guess_matrix,fine_g='tanh', mu=0.01, stabilization=True, sample_size=1, fine_tanh=1, max_it=500000, max_it_fine=1000, failures=5000, limit=0.0005, verbose=True, whitened=True)

#self.source_list=mdp.fastica(self.source_list, approach='symm',g='gaus', guess=guess_matrix,fine_g='gaus', mu=0.1, stabilization=0.01, sample_size=1, fine_tanh=1, max_it=5000000, max_it_fine=100, failures=500, limit=0.001, verbose=False, whitened=False, white_comp=self.nout, white_parm={'svd':True})
        #ficad=mdp.nodes.FastICANode(approach='defl',g='pow3', guess=guess_matrix,fine_g='pow3', mu=0.1, stabilization=True, sample_size=1, fine_tanh=1, max_it=1000000, max_it_fine=100000, failures=5000, limit=100000, verbose=True, whitened=True)
        ficad.train(whitened_source)
        
        self.source_list=ficad.execute(whitened_source)
        f=open("ica_matrix.dat","w")
        
        #ficad.get_projmatrix().tofile(open("ica_matrix.dat","w"),sep=" ")
        ica_matrix=ficad.get_projmatrix()
        for i in range(0,ica_matrix.shape[0]):
            for j in range(0,ica_matrix.shape[1]):
                f.write("%s "%ica_matrix[i,j])
            f.write("\n")
        f.close()

        print ica_matrix

        #        proj_matrix=mdp.utils.mult(whiten_matrix,ica_matrix)
        #        proj_matrix.tofile("proj_matrix.dat",sep=" ")
        
        #        rec_matrix=mdp.utils.inv(proj_matrix)
        #print rec_matrix

#        for i in range(0,rec_matrix.shape[0]):
#            indices.append(rec_matrix[i,0]/rec_matrix[i,-1]);
#        indices=numpy.abs(indices)
#        indices.sort()
#        print indices
        result=[]
#        print self.source_list.shape
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

    #print sys.argv[1]==none
    if sys.argv[1]=="none":
        mask=None
    else:
        mask=pyfits.open(sys.argv[1])[0].data
        
    
    
    nout=int(sys.argv[2])
    out_prefix=sys.argv[3+nsources]
#    print "out prefix:",out_prefix
    mi=mask_icaer(nout)
    mi.mask=(mask.astype('int'))
    for i in range(3,3+nsources):
#        print "loading:",sys.argv[i]
        img=pyfits.open(sys.argv[i])[0].data
        if mask==None:
            mask=numpy.ones(img.shape)
        mi.append_source(img)

#    print "Performing ICA"
    
    sources=mi.perform()

    for i in range(0,nout):
        outname=out_prefix+`i+1`+".fits"
#        print "Writing ",outname
        pyfits.HDUList([pyfits.PrimaryHDU(sources[i])]).writeto(outname,clobber=True)


        
