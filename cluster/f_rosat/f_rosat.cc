#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cmath>
//#include "calc_distance.h"
//#include "ddivid.h"
using namespace std;


static double fconv_fact[][4]=
  {
#include "fconv_table.inc"
  };

static int n_z;
static int n_T;
static int n_row=sizeof(fconv_fact)/(4*sizeof(double));

static vector<double> T_vec;
static vector<double> z_vec;

double current_T;

static class initlizer
{
public:
  initlizer()
  {
    n_z=1;
    n_T=1;
    double z=fconv_fact[0][0];
    double T=fconv_fact[0][1];
    z_vec.push_back(z);
    T_vec.push_back(T);
    for(int i=0;i<n_row;++i)
      {
	if(z<fconv_fact[i][0])
	  {
	    ++n_z;
	    z=fconv_fact[i][0];
	    z_vec.push_back(z);
	  }
	if(T<fconv_fact[i][1])
	  {
	    ++n_T;
	    T=fconv_fact[i][1];
	    T_vec.push_back(T);
	  }
      }
  }
}_initlizer;

static int grn(int z_index,int T_index)
{
  return z_index*n_T+T_index;
}

static int get_T_index(double T)
{
  assert(T>=T_vec[0]);
  for(int i=0;i<n_T-1;++i)
    {
      if(T>=T_vec[i]&&T<T_vec[i+1])
	return i;
    }
  return -1;
}

static int get_z_index(double z)
{
  assert(z>=z_vec[0]);
  for(int i=0;i<n_z-1;++i)
    {
      if(z>=z_vec[i]&&z<z_vec[i+1])
	return i;
    }
  return -1;
}

static double get_bolo_flux(int iz,int iT)
{
  return fconv_fact[grn(iz,iT)][2];
}

static double get_soft_flux(int iz,int iT)
{
  return fconv_fact[grn(iz,iT)][3];
}

double get_bolo_flux(double z,double T)
{
  int iT=get_T_index(T);
  int iz=get_z_index(z);
  double f00=get_bolo_flux(iz,iT);
  double f01=get_bolo_flux(iz,iT+1);
  double f10=get_bolo_flux(iz+1,iT);
  
  double T0=T_vec[iT];
  double T1=T_vec[iT+1];
  double z0=z_vec[iz];
  double z1=z_vec[iz+1];

  return f00+
    (f01-f00)/(T1-T0)*(T-T0)+
    (f10-f00)/(z1-z0)*(z-z0);
   
}

double get_soft_flux(double z,double T)
{
  int iT=get_T_index(T);
  int iz=get_z_index(z);
  double f00=get_soft_flux(iz,iT);
  double f01=get_soft_flux(iz,iT+1);
  double f10=get_soft_flux(iz+1,iT);
  
  double T0=T_vec[iT];
  double T1=T_vec[iT+1];
  double z0=z_vec[iz];
  double z1=z_vec[iz+1];

  return f00+
    (f01-f00)/(T1-T0)*(T-T0)+
    (f10-f00)/(z1-z0)*(z-z0);
   
}
