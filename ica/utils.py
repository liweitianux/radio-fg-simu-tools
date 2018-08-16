#!/usr/bin/env python

import numpy


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
    #        print mx.shape
    #        print size1
    for i in range(0,mask.shape[0]):
        for j in range(0,mask.shape[1]):
            if mask[i,j]!=0:
                result[i,j]=mx[n]
                n+=1
    return result
