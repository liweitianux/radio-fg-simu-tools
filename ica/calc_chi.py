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


sbchi_file=open('sb_chi.qdp','w')

sbchi_file.write('read serr 2\n')
sbchi_file.write('skip single\n')
sbchi_file.write('ma 3 on 1\n')
sbchi_file.write('li on 2\n')
sbchi_file.write('la x Radius (arcmin)\n')
sbchi_file.write('la y chi\n')
sbchi_file.write('la f\n')
sbchi_file.write('r y -2 2\n')
sbchi_file.write('time off\n')

chi_sq=0

for i in range(0,len(x_list)):
    sbchi_file.write('%s %s %s\n'%(x_list[i],(y_list[i]-ymodel_list[i])/yerr_list[i],1))
    chi_sq+=((y_list[i]-ymodel_list[i])/yerr_list[i])**2
    
chi_sq/=len(x_list)
print "sb chi^2=",chi_sq


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


specchi_file=open('spec_total_chi.qdp','w')

specchi_file.write('read serr 2\n')
specchi_file.write('skip single\n')
specchi_file.write('ma 3 on 1\n')
specchi_file.write('li on 2\n')
specchi_file.write('la x Radius (arcmin)\n')
specchi_file.write('la y chi\n')
specchi_file.write('la f\n')
specchi_file.write('r y -2 2\n')
specchi_file.write('time off\n')

chi=0
for i in range(0,len(x_list)):
    specchi_file.write('%s %s %s\n'%(x_list[i],(y_list[i]-ymodel_list[i])/yerr_list[i],1))
    chi_sq+=((y_list[i]-ymodel_list[i])/yerr_list[i])**2

chi_sq/=len(x_list)
print "spec_total chi^2=",chi_sq

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


specchi_file=open('spec_peak_chi.qdp','w')

specchi_file.write('read serr 2\n')
specchi_file.write('skip single\n')
specchi_file.write('ma 3 on 1\n')
specchi_file.write('li on 2\n')
specchi_file.write('la x Radius (arcmin)\n')
specchi_file.write('la y chi\n')
specchi_file.write('la f\n')
specchi_file.write('r y -2 2\n')
specchi_file.write('time off\n')

for i in range(0,len(x_list)):
    specchi_file.write('%s %s %s\n'%(x_list[i],(y_list[i]-ymodel_list[i])/yerr_list[i],1))
