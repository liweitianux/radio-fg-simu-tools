#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <cassert>
#include <vector>
#include "adapt_trapezoid.h"
#include "coscalc.h"
//calc_distance
//usage:
//calc_distance z

using namespace std;

const double cm=1;
const double s=1;
const double km=1000*100;
const double Mpc=3.08568e+24*cm;
const double kpc=3.08568e+21*cm;
const double yr=365.*24.*3600.;
const double Gyr=1e9*yr;
double H=70.*km/s/Mpc;
const double c=299792458.*100.*cm;
double omega_m=0.3;
double omega_l=0.7;
double omega_k=1-omega_m-omega_l;
const double arcsec2arc_ratio=1./60/60/180*3.1415926;

double E1(double z);

static class int_E_table
{
public:
  const double z_max;
  const double z_min;
  const double delta_z;
  vector<double> E_table;
public:
  int_E_table()
    :z_max(22.01),z_min(.01),delta_z(.001)
  {
    E_table.resize((z_max-z_min)/delta_z);
    for(size_t i=0;i<E_table.size();++i)
      {
	double z=z_min+i*delta_z;
	E_table[i]=adapt_trapezoid(E1,0.,z,1E-3);
      }
    cerr<<"E table initialized"<<endl;
  }
}E_table;


double calc_int_E(double z)
{
  if(z>E_table.z_min&&z<E_table.z_max)
    {
      return E_table.E_table[(z-E_table.z_min)/E_table.delta_z];
    }
  return adapt_trapezoid(E1,0.,z,1E-3);
}


void update_const(double H0,double Om,double Ol)
{
  H=H0;
  omega_m=Om;
  omega_l=Ol;  
  omega_k=1-omega_m-omega_l;
}

double E(double z)
{
  return sqrt(omega_m*(1+z)*(1+z)*(1+z)+omega_k*(1+z)*(1+z)+omega_l);
  //return 1/sqrt(1+z*omega_m+(1/(1+z)/(1+z)-1)*omega_l);
}


double E1(double z)
{
  return 1/E(z);
}

double calc_comoving_distance(double z)
{
  //return c/H*adapt_trapezoid(E1,0.,z,1E-1);
  return c/H*calc_int_E(z);
}


double calc_angular_distance(double z)
{
  return calc_comoving_distance(z)/(1+z);
}

double calc_dvc_dz(double z,double dO)
{
  double da=calc_angular_distance(z);
  return c/H*(1+z)*(1+z)*da*da/E(z)*dO;
}
