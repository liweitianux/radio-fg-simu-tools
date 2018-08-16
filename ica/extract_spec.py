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
    print i
    profile_sum.append(numpy.zeros(rmax+1))
    profile_cnt.append(numpy.zeros(rmax+1))
#    num_list.append(i)



c=3e10

kb=1.38E-16

pixel_area=(10/180.*math.pi/1024)**2


def Tb2flux(Tb,nu_MHz):
    nu=nu_MHz*1E6
    return kb*Tb*2*nu**2/c**2*pixel_area

fout=open('spec.qdp','w')

for k in range(0,len(freq_list)):
    fname=fprefix+`freq_list[k]`+'.fits'
    img=pyfits.open(fname)[0].data
    center_ij=img_utils.max_element(img)
    print center_ij

    for i in range(-rmax,rmax+1):
        for j in range(-rmax,rmax+1):
            r=int(math.sqrt(i**2+j**2))
            if img[center_ij[0]+i,center_ij[1]+j] >0 and r<rmax:
                profile_sum[k][int(r)]+=img[center_ij[0]+i,center_ij[1]+j]
                profile_cnt[k][int(r)]+=1

    for i in range(0,len(profile_sum[k])):
        if profile_cnt[k][i]>0:
            profile_sum[k][i]/=profile_cnt[k][i]
            if profile_sum[k][i]<profile_sum[k][0]/10.:
                rmax=min(rmax,i)
#            print i,profile_sum[i]

    
print "i j r=\n"
print center_ij[0],center_ij[1],rmax

rmax=1

errors=open('error.dat').readlines()

for k in range(0,len(freq_list)):
    total_flux=0
#    print profile_sum[k][0]
    for i in range(0,rmax):
        if profile_cnt[k][i]>0:
 #           print profile_sum[i]
            total_flux+=2*math.pi*(i+0.5)*Tb2flux(profile_sum[k][i],freq_list[k])/1E-23
            
            #total_flux=profile_sum[k][0]
    #            print profile_sum[i],Tb2flux(Tb,freq_list[k])/1E-23
    fout.write("%E %E %E\n"%(freq_list[k],total_flux,total_flux*float(errors[k])))

#    print 'no no no'
fout.write("no no no\n")
