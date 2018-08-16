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


rmax=0

if len(sys.argv)!=5:
    print "Usage:",sys.argv[0]," <input fits> <output qdp> <rmax> <error ratio>"
    sys.exit(-1)


rmax=float(sys.argv[3])
fout=open(sys.argv[2],'w')

img=pyfits.open(sys.argv[1])[0].data



peak_i=0
peak_j=0
max_value=0

for i in range(0,img.shape[0]):
    for j in range(0,img.shape[1]):
        if max_value<img[i,j]:
            max_value=img[i,j]
            peak_i=i
            peak_j=j



profile_sum=numpy.zeros([rmax])
profile_cnt=numpy.zeros([rmax])

for di in range(-rmax,rmax+1):
    for dj in range(-rmax,rmax+1):
        if img[di+peak_i,dj+peak_j]!=0:
            r=math.sqrt(di**2+dj**2)
            r=int(r)
            if r<rmax:
                profile_sum[r]+=img[di+peak_i,dj+peak_j]
                profile_cnt[r]+=1

error_ratio=float(sys.argv[4])




for i in range(0,rmax):
    if profile_cnt[i]>0:
        profile_sum[i]/=profile_cnt[i]

fout.write('read serr 2\n')
fout.write('skip single\n')
fout.write('ma 3 on 1\n')
fout.write("r y %s %s\n"%(min(profile_sum)*(1-error_ratio),max(profile_sum)*(1+error_ratio)))



for i in range(0,rmax):
    fout.write("%s %s %s\n"%(i*0.2929,profile_sum[i],error_ratio*profile_sum[i]))

