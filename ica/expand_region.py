#!/usr/bin/env python

import re
import numpy
import scipy
import mdp
import sys
import math
import exceptions
import img_utils
import pyfits


if len(sys.argv)!=3:
    print "Usage:", sys.argv[0]," <input region file> <output prefix>"
    sys.exit(-1)




###########
delta_r=10
delta_angle=5
###########

cnt=1
for i in open(sys.argv[1]):
    if re.match('circle',i):

        splited=re.findall(r"(\S+)\((\S+),(\S+),(\S+)\)",i.strip())

        if len(splited)>0:
            divname=sys.argv[2]+'div_'+`cnt`+'.reg'
            expname=sys.argv[2]+'exp_'+`cnt`+'.reg'
            
            divfile=open(divname,'w')
            expfile=open(expname,'w')
            
            
            jj=float(splited[0][1])
            ii=float(splited[0][2])
            r=float(splited[0][3])
            r+=delta_r;
            r=min(r,512)
            divfile.write(i+'\n')
            expfile.write("circle(%s,%s,%s)\n"%(str(jj),str(ii),str(r)))
            cnt+=1

    elif re.match('pie',i):
        splited=re.findall(r"(\S+)\((\S+),(\S+),(\S+),(\S+),(\S+),(\S+)\)",i.strip())
        if len(splited)>0:
            divname=sys.argv[2]+'div_'+`cnt`+'.reg'
            expname=sys.argv[2]+'exp_'+`cnt`+'.reg'
            
            divfile=open(divname,'w')
            expfile=open(expname,'w')
            
            
            jj=float(splited[0][1])
            ii=float(splited[0][2])
            ri=float(splited[0][3])
            ro=float(splited[0][4])
            a1=float(splited[0][5])
            a2=float(splited[0][6])

            ri-=delta_r
            ro+=delta_r
            
            ri=max(0,ri)
            ro=min(ro,512)
            a1-=delta_r/ri*180./math.pi
            a2+=delta_r/ri*180./math.pi
            divfile.write(i+'\n')
            expfile.write("pie(%s,%s,%s,%s,%s,%s)\n"%(str(jj),str(ii),str(ri),str(ro),str(a1),str(a2)))
            cnt+=1
