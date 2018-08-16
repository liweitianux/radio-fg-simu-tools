#!/usr/bin/env python


import numpy
import scipy
import mdp
import sys
import math
import exceptions
import pyfits
import img_utils

#num_list=[]
freq_list=[65,80,95,110,125,140,155,170,185]


if len(sys.argv)!=3:
    print "Usage:",sys.argv[0]," <rmax> <prefix>"
    sys.exit(-1)

rmax=float(sys.argv[1])
fprefix=sys.argv[2]

profile_sum=[]
profile_cnt=[]

for i in range(0,len(freq_list)):
    profile_sum.append(numpy.zeros(rmax+1))
    profile_cnt.append(numpy.zeros(rmax+1))
#    num_list.append(i)



c=3e10

kb=1.38E-16

pixel_area=(10/180.*math.pi/2048)**2


def Tb2flux(Tb,nu_MHz):
    nu=nu_MHz*1E6
    return kb*Tb*2*nu**2/c**2*pixel_area


errors=open('error.dat').readlines()

for k in range(0,len(freq_list)):
    fname=fprefix+`freq_list[k]`+'.fits'
    img=pyfits.open(fname)[0].data
    center_ij=img_utils.max_element(img)
    print freq_list[k],img[center_ij],float(errors[k])*img[center_ij]

print "no no no"
