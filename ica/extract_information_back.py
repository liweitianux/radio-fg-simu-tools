#!/usr/bin/env python
from img_utils import *
import pyfits
#extracting spectrum

import math
c=3e10
kb=1.38E-16
pixel_area=(10/180.*math.pi/1024)**2

def pixel2arcmin(x):
    return x*0.5859


def Tb2flux(Tb,nu_MHz):
    nu=nu_MHz*1E6
    return kb*Tb*2*nu**2/c**2*pixel_area


mask=pyfits.open('mask2.fits')[0].data
img65=pyfits.open('65.fits')[0].data
img_model=pyfits.open('model_65.fits')[0].data

peak_position=max_element(img_model)

peak_value=img65[peak_position]


rmax=0
for r in range(1,200):
    circle_pixels=draw_circle(peak_position[0],peak_position[1],r)
    sum=0
    for i in circle_pixels:
        sum+=img65[i]
    sb=sum/len(circle_pixels)
    if sb<=peak_value/10:
        rmax=r
        break
    print r,sb
print r

#extracting 65 MHz SB

sb_file=open('sb.qdp','w')
sb_file.write('read serr 2\n')
sb_file.write('skip single\n')
sb_file.write('ma 3 on 1\n')
sb_file.write('li on 2\n')
sb_file.write('la x Radius (arcmin)\n')
sb_file.write('la y Brightness temperature (K)\n')
sb_file.write('la f\n')
sb_file.write('time off\n')

sberror=0

for i in range(1,rmax):
    circle_pixels=draw_circle(peak_position[0],peak_position[1],i)
    sum=0
    sb2sum=0
    max_sb=0
    min_sb=1e99
    cnt=0
    for j in circle_pixels:
        if mask[j] >0:
            sum+=img65[j]
            max_sb=max(max_sb,img65[j])
            min_sb=min(min_sb,img65[j])
            sb2sum+=img65[j]**2
            cnt+=1
    sb=sum/cnt
    #sberror+=(max_sb-min_sb)/(max_sb+min_sb)

    sbsigma=math.sqrt(sb2sum/cnt-sb**2)
    sberror+=sbsigma/(max_sb+min_sb)    
    sb_file.write(`pixel2arcmin(i)`+" "+`sb`+" "+`sbsigma`+'\n')

sberror/=rmax
    
sb_file.write('no no no\n')


for i in range(1,rmax):
    circle_pixels=draw_circle(peak_position[0],peak_position[1],i)
    sum=0
    max_sb=0
    min_sb=1e99
    for j in circle_pixels:
        sum+=img_model[j]
        max_sb=max(max_sb,img65[j])
        min_sb=min(min_sb,img65[j])
    sb=sum/len(circle_pixels)
    sb_file.write(`pixel2arcmin(i)`+" "+`sb`+" "+`0`+'\n')

#extracting the spectrum

sb_file.close()
sberrs={}

for f in [65,80,95,110,125,140,155,170,185]:
    img=pyfits.open('%s.fits'%(`f`))[0].data
    sberror1=0
    for i in range(1,rmax):
        circle_pixels=draw_circle(peak_position[0],peak_position[1],i)
        sum=0
        sb2sum=0
        max_sb=0
        min_sb=1e99
        cnt=0
        for j in circle_pixels:
            if mask[j]>0:
                sum+=img65[j]
                max_sb=max(max_sb,img[j])
                min_sb=min(min_sb,img[j])
                sb2sum+=img65[j]**2
                cnt+=1
        sb=sum/cnt
    #sberror+=(max_sb-min_sb)/(max_sb+min_sb)

        sbsigma=math.sqrt(sb2sum/cnt-sb**2)
        sberror1+=sbsigma/(max_sb+min_sb)    
    sberror1/=rmax
    print sberror1
    sberrs[int(f)]=sberror1
       


#####
spec_total_file=open('spec_total.qdp','w')

spec_total_file.write('read serr 2\nskip single\n')
spec_total_file.write('ma 2 on 1\nli on 2\n')
spec_total_file.write('r x 60 200\n')
spec_total_file.write('log\n')
spec_total_file.write('la x Frequency (MHz)\n')
spec_total_file.write('la y Flux (Jy)\n')
spec_total_file.write('la f\n')
spec_total_file.write('time off\n')
for freq in range(65,186,15):
    print freq
    img=pyfits.open(`freq`+'.fits')[0].data
    sum=0
    for i in range(peak_position[0]-rmax,peak_position[0]+rmax):
        for j in range(peak_position[1]-rmax,peak_position[1]+rmax):
            r=math.sqrt((i-peak_position[0])**2+(j-peak_position[1])**2)
            
            if r<rmax:
                sum+=img[i,j]
    spec_total_file.write(`freq`+' '+`Tb2flux(sum,freq)/1e-23`+' '+`Tb2flux(sum,freq)*sberrs[freq]/1e-23`+'\n')
    #spec_total_file.write(`freq`+' '+`sum`+' '+`sum*sberror`+'\n')

                
                
spec_total_file.write('no no no\n')

for freq in [65,80,95,110,125,140,155,170,185]:
    img=pyfits.open('model_'+`freq`+'.fits')[0].data
    sum=0
    for i in range(peak_position[0]-rmax,peak_position[0]+rmax):
        for j in range(peak_position[1]-rmax,peak_position[1]+rmax):
            r=math.sqrt((i-peak_position[0])**2+(j-peak_position[1])**2)
            
            if r<rmax:
                sum+=img[i,j]
    spec_total_file.write(`freq`+' '+`Tb2flux(sum,freq)/1e-23`+' '+'0'+'\n')
    #spec_total_file.write(`freq`+' '+`sum`+' '+'0'+'\n')

spec_total_file.close()


spec_peak_file=open('spec_peak.qdp','w')
spec_peak_file.write('read serr 2\nskip single\n')
spec_peak_file.write('ma 2 on 1\nli on 2\n')
spec_peak_file.write('r x 60 200\n')
spec_peak_file.write('log\n')
spec_peak_file.write('la x Frequency (MHz)\n')
spec_peak_file.write('la y Brightness Temperature (K)\n')
spec_peak_file.write('la f\n')
spec_peak_file.write('time off\n')

for freq in range(65,186,15):
    img=pyfits.open(`freq`+'.fits')[0].data
    kT_peak=img[peak_position]
    spec_peak_file.write(`freq`+' '+`kT_peak`+' '+`sberrs[freq]*kT_peak`+'\n')


spec_peak_file.write('no no no\n')
for freq in range(65,186,15):
    img=pyfits.open('model_'+`freq`+'.fits')[0].data
    kT_peak=img[peak_position]
    spec_peak_file.write(`freq`+' '+`kT_peak`+' '+`0`+'\n')


