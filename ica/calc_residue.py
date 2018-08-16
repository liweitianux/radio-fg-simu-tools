#!/usr/bin/env python

sb_file=open('sb.qdp')

flag=0

x_list=[]
y_list=[]
yerr_list=[]
ymodel_list=[]

for i in sb_file:
    if i.strip()=='time off':
        flag=1
        continue
    if i.strip()=='no no no':
        flag=2
        continue
    if flag==1:
        dummy=0

        x,y,yerr=i.split()
        x_list.append(float(x))
        y_list.append(float(y))
        yerr_list.append(float(yerr))
        
    if flag==2:
        x,y,yerr=i.split()
        ymodel_list.append(float(y))


sbresidue_file=open('sb_residue.qdp','w')

sbresidue_file.write('read serr 2\n')
sbresidue_file.write('skip single\n')
sbresidue_file.write('ma 3 on 1\n')
sbresidue_file.write('li on 2\n')
sbresidue_file.write('la x Radius (arcmin)\n')
sbresidue_file.write('la y residue (K)\n')
sbresidue_file.write('la f\n')
sbresidue_file.write('time off\n')

for i in range(0,len(x_list)):
    sbresidue_file.write('%s %s %s\n'%(x_list[i],(y_list[i]-ymodel_list[i]),0))


################
spec_file=open('spec_total.qdp')

flag=0

x_list=[]
y_list=[]
yerr_list=[]
ymodel_list=[]

for i in spec_file:
    if i.strip()=='time off':
        flag=1
        continue
    if i.strip()=='no no no':
        flag=2
        continue
    if flag==1:
        dummy=0

        x,y,yerr=i.split()
        x_list.append(float(x))
        y_list.append(float(y))
        yerr_list.append(float(yerr))
        
    if flag==2:
        x,y,yerr=i.split()
        ymodel_list.append(float(y))


specresidue_file=open('spec_total_residue.qdp','w')

specresidue_file.write('read serr 2\n')
specresidue_file.write('skip single\n')
specresidue_file.write('ma 3 on 1\n')
specresidue_file.write('li on 2\n')
specresidue_file.write('la x Radius (arcmin)\n')
specresidue_file.write('la y residue (Jy)\n')
specresidue_file.write('la f\n')
specresidue_file.write('time off\n')

for i in range(0,len(x_list)):
    specresidue_file.write('%s %s %s\n'%(x_list[i],(y_list[i]-ymodel_list[i]),0))


################
spec_file=open('spec_peak.qdp')

flag=0

x_list=[]
y_list=[]
yerr_list=[]
ymodel_list=[]

for i in spec_file:
    if i.strip()=='time off':
        flag=1
        continue
    if i.strip()=='no no no':
        flag=2
        continue
    if flag==1:
        dummy=0

        x,y,yerr=i.split()
        x_list.append(float(x))
        y_list.append(float(y))
        yerr_list.append(float(yerr))
        
    if flag==2:
        x,y,yerr=i.split()
        ymodel_list.append(float(y))


specresidue_file=open('spec_peak_residue.qdp','w')

specresidue_file.write('read serr 2\n')
specresidue_file.write('skip single\n')
specresidue_file.write('ma 3 on 1\n')
specresidue_file.write('li on 2\n')
specresidue_file.write('la x Radius (arcmin)\n')
specresidue_file.write('la y residue (K)\n')
specresidue_file.write('la f\n')
specresidue_file.write('time off\n')

for i in range(0,len(x_list)):
    specresidue_file.write('%s %s %s\n'%(x_list[i],(y_list[i]-ymodel_list[i]),0))
