#include <cmath>
#include "coscalc.hpp"
#include "adapt_trapezoid.h"
#include <iostream>
#include <cassert>
using namespace std;

double cm=1;
double gram=1;
double sec=1;


double erg=gram*cm*cm/sec/sec;
double meter=100*cm;
double kg=gram*1000;
double km=meter*1000;


double hr=3600*sec;
double day=24*hr;
double yr=365.2422*day;
double kpc=3.08568e+21*cm;
double Mpc=1000*kpc;
double J=1e7*erg;


double mp=1.6726215813e-27*kg;
double me=9.1093818872e-31*kg;
double eV=1.602176462e-19*J;
double keV=1000*eV;



double H0=71*km/sec/Mpc;
double c=2.99792458E8*meter/sec;


void init_unit_system(int us)
{
  if(us==CGS)
    {
      cm=1;
      sec=1;
      gram=1;
    }
  else if(us==SI)
    {
      cm=.01;
      sec=1;
      gram=1E-3;
    }
  else
    {
      cerr<<"unit system undefined!"<<endl;
      assert(0);
    }
  double erg=gram*cm*cm/sec/sec;
  double meter=100*cm;
  double kg=gram*1000;
  double km=meter*1000;
  
  
  double hr=3600*sec;
  double day=24*hr;
  double yr=365.2422*day;
  double kpc=3.08568e+21*cm;
  double Mpc=1000*kpc;
  double J=1e7*erg;
  
  
  double mp=1.6726215813e-27*kg;
  double me=9.1093818872e-31*kg;
  double eV=1.602176462e-19*J;
  double keV=1000*eV;
  
  
  
  double H0=71*km/sec/Mpc;
  double c=2.99792458E8*meter/sec;
  
  
}


double E(double z,double omega_m,double omega_l)
{
  double omega_k=1-omega_m-omega_l;
  return sqrt(omega_m*(1+z)*(1+z)*(1+z)+omega_k*(1+z)*(1+z)+omega_l);
}

static double int_E(double z,void* param)
{
  double omega_m=*((double*)param+0);
  double omega_l=*((double*)param+1);

  return 1/E(z,omega_m,omega_l);
}


double calc_dc(double z,double omega_m,double omega_l)
{
  double p[2];
  p[0]=omega_m;
  p[1]=omega_l;
  return c/H0*adapt_trapezoid(int_E,0.,z,1E-4,p);
}

double calc_dm(double z,double omega_m,double omega_l)
{
  double omega_k=1-omega_m-omega_l;
  double Dh=c/H0;
  double Dc=calc_dc(z,omega_m,omega_l);
  if(omega_k>0)
    {
      return Dh/sqrt(omega_k)*sinh(sqrt(omega_k)*Dc/Dh);
    }
  else if(omega_k<0)
    {
      return Dh/sqrt(omega_k)*sin(sqrt(omega_k)*Dc/Dh);
    }
  else
    {
      return Dc;
    }
  
}

double calc_da(double z,double omega_m,double omega_l)
{
  return calc_dm(z,omega_m,omega_l)/(1+z);
}

double calc_dl(double z,double omega_m,double omega_l)
{
  return calc_dm(z,omega_m,omega_l)*(1+z);
}
