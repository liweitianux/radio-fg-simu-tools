import scipy
import scipy.special
import scipy.integrate
import math
from math import pi

e=4.8E-10
me=9.1E-28
c=2.99792458E10

gmin=1
gmax=1E6

"""
Package used to calculate the decay time of synchrotron emission.
"""


def kv53(x):
    return scipy.special.kv(5/3.,x)

Fn=.007*scipy.integrate.quad(kv53,.007,1000)[0]

def F(x):
    """F(v/vc)"""
    if x<10 and x>0.007:
        return x*scipy.integrate.quad(kv53,x,1000)[0]
    elif x<=.007:
        return Fn*(x/.007)**(1/3.)
    elif x>=10:
        return (pi/2*x)**.5*math.exp(-x)
    return 0

def nu_L(B):
    """Calculate nu_L based on magnetic field strength"""
    return e*B/(2*pi*me*c)

def nu_c(B,gm):
    """Calculate the critical frequency"""
    return 1.5*nu_L(B)*gm**2

def single_spec(B,gm,nu):
    """The spectrum of one single electron"""
    x=nu/nu_c(B,gm)
    #    print x
    return 2*pi*3**.5*e**2*nu_L(B)/c*F(x)

def N(b,B,gm,n,N0,t):
    """Time evolution of a group of electrons with power law distribution
    at the original time. 
    Parameters:
    b: The decay rate, calculated with 3E-8/(8 pi)*B^2
    B: Local magnetic field
    gm: gamma
    n: The initial power law index
    N0: Normalization coefficient
    t: time in unit of sec
    """
    if 1<=b*gm*t:
        return 0
    else:
        return N0*gm**(-n)*(1-b*gm*t)**(n-2)
    
def int_func(gm,*args):
    b=args[0]
    B=args[1]
    n=args[2]
    nu=args[3]
    N0=args[4]
    t=args[5]
    return N(b,B,gm,n,N0,t)*single_spec(B,gm,nu)

points=[gmin,gmax/1E6,gmax/1E5,gmax/1E4,gmax/1E3,gmax/1E2,gmax]

def int_spec(b,B,n,nu,N0,t):
    """
    The integrated spectrum of a group of electrons at time t
    Parameters: 
    b: The decay rate, calculated with 3E-8/(8 pi)*B^2
    B: Local magnetic field
    n: The initial power law index
    nu: Radio frequency in unit Hz
    N0: The normalization coefficient of electron distribution
    t: time in unit second
    """
    epsabs=1E-35*N0
    result=scipy.integrate.quad(int_func,gmin,gmax,args=(b,B,n,nu,N0,t),epsabs=epsabs,points=points);
    while result[1]>result[0]*.01:
        epsabs*=1E-2
        print epsabs
        result=scipy.integrate.quad(int_func,gmin,gmax,args=(b,B,n,nu,N0,t),epsabs=epsabs,points=points);
        
   
    return result[0]        

def calc_b(B1,B2):
    """
    The help function to calculate the decay rate from cosmological B and local B
    """ 
    return 3E-8/8/pi*(B1**2+B2**2)

def plot_F(t_list=[0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0],v1=5E6,v2=2E9):
    """
    Plot a series of spectra at different time
    Parameters:
    t_list: a sequence of time in unit Gyr
    v1: The lower limit of frequency
    v2: The upper limit of frequency
    """
    import pylab
    for t in t_list:
        x=[]
        y=[]
        i=v1
        while i<v2:
            yt=int_spec(calc_b(0.2E-6,3.2E-6*(1+0.)*(1+0.)),0.2E-6,3,i,1,t*1E9*365*24*3600)
            if yt<=0:
                break
            x.append(i)
            y.append(yt)
            print i
            i*=1.2
        pylab.loglog(x,y)

    
            
def dump_qdp(outfile=None,t_list=[0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0],v1=5E6,v2=2E9):
    """
    Plot a series of spectra at different time
    Parameters:
    t_list: a sequence of time in unit Gyr
    v1: The lower limit of frequency
    v2: The upper limit of frequency
    """
    if outfile==None:
        outfile=open("dump.qdp",'w')

    outfile.write("skip single\n")
    outfile.write("log\n")

    for t in t_list:
        x=[]
        y=[]
        i=v1
        while i<v2:
            yt=int_spec(calc_b(0.2E-6,3.2E-6*(1+0.)*(1+0.)),0.2E-6,3,i,1,t*1E9*365*24*3600)
            if yt<=0:
                break
#            x.append(i)
#            y.append(yt)
            outfile.write("%s %s\n"%(str(i),str(yt*1E40)))
            print i
            i*=1.2
        outfile.write("no no no\n")
        outfile.flush()
    outfile.close()
