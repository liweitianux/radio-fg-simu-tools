#!/usr/bin/env python


import numpy
import scipy
import mdp
import sys
import math
import exceptions
import pyfits
import img_utils

num_list=[]
freq_list=[65,80,95,110,125,140,155,170,185]

rmax=0


fprefix="/data2/simu_diffuse/cluster/Timg_cluAB_"
fsuffix=".fits"


if len(sys.argv)!=4:
    print "Usage:",sys.argv[0]," i j r"

    sys.exit(-1)


center_i=int(sys.argv[1])
center_j=int(sys.argv[2])

rmax=int(sys.argv[3])

profile_sum=[]
profile_cnt=[]



c=3e10

kb=1.38E-16

pixel_area=(10/180.*math.pi/2048)**2


def Tb2flux(Tb,nu_MHz):
    nu=nu_MHz*1E6
    return kb*Tb*2*nu**2/c**2*pixel_area


for n in freq_list:
    img=pyfits.open(fprefix+`n`+fsuffix)[0].data
    total_flux=0
    for i in range(-rmax,rmax+1):
        for j in range(-rmax,rmax+1):
            if i*i+j*j<=rmax**2:
                total_flux+=Tb2flux(img[center_i+i,center_j+j],n)/1e-23
                #print Tb2flux(img[center_i+i,center_j+j],n)/1e-23
    #total_flux=img[center_i,center_j]
    print n,total_flux

