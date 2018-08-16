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


if len(sys.argv)!=4:
    print "Usage:", sys.argv[0]," <region file> <img> <output mask>"
    sys.exit(-1)

mask=numpy.zeros(pyfits.open(sys.argv[2])[0].data.shape,'int')
for i in open(sys.argv[1]):
    if re.match('circle',i):
        splited=re.findall(r"(\S+)\((\S+),(\S+),(\S+)\)",i.strip())

        if len(splited)>0:
            jj=float(splited[0][1])
            ii=float(splited[0][2])
            r=float(splited[0][3])
#            print ii,jj,r
            for n in img_utils.make_round_region(mask,ii,jj,r):
                mask[n]=1
    elif re.match('pie',i):
        splited=re.findall(r"(\S+)\((\S+),(\S+),(\S+),(\S+),(\S+),(\S+)\)",i.strip())
        if len(splited)>0:
            jj=float(splited[0][1])
            ii=float(splited[0][2])
            ri=float(splited[0][3])
            ro=float(splited[0][4])
            a1=float(splited[0][5])
            a2=float(splited[0][6])

#            print ii,jj,ri,ro,a1,a2
            for n in img_utils.make_pie_region(mask,ii,jj,ri,ro,a1,a2):
                mask[n]=1
    elif re.match('box',i):
        splited=re.findall(r"(\S+)\((\S+),(\S+),(\S+),(\S+),(\S+)\)",i.strip())
        if len(splited)>0:
            jj=float(splited[0][1])
            ii=float(splited[0][2])
            w=float(splited[0][3])
            h=float(splited[0][4])
#            print ii,jj,w,h
            for n in img_utils.make_box_region(mask,ii,jj,h,w):
                mask[n]=1
pyfits.HDUList(pyfits.PrimaryHDU(mask)).writeto(sys.argv[3],clobber=True)


