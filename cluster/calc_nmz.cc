#include <fstream>
#include <iostream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <cassert>
#include "adapt_trapezoid.h"
#define IMASS_FUNC_A

using namespace std;

double lower_mass;
const double upper_redshift=3.0;
double cm=1;
double s=1;
double km=1000*100;
double Mpc=3.08568e+24;
double kpc=3.08568e+21;
double yr=365.*24.*3600.;
double Gyr=1e9*yr;
double H=71.*km/s/Mpc;//Hubble constant changed 70->71 2009/08/20
const double c=299792458.*100;
const double omega_m=0.27;//constant changed 2008/08/20
const double omega_l=0.73;
const double arcsec2arc_ratio=1./60/60/180*3.1415926;
const int number_z=149;
double da;

double E(double z)
{
  double omega_k=1-omega_m-omega_l;
  return sqrt(omega_m*(1+z)*(1+z)*(1+z)+omega_k*(1+z)*(1+z)+omega_l);
}

double f_dist(double z)
{
  //  return 1/((1+z)*sqrt(1+omega_m*z-omega_l*(1-1/((1+z)*(1+z)))));
  return 1/E(z);
}

double calc_H(double z)
{
  //  return H/f_dist(z);
  return H*E(z);
}

double calc_angular_distance(double z)
{
  return c/H*adapt_trapezoid(f_dist,0.,z,1e-4)/(1+z);
}

int main(int argc,char* argv[])
{ 
  if (argc != 3)
    {
        cerr << "Usage:" << endl;
        cerr << "    " << argv[0] << "LowerMass dndmdata" << endl;
        exit(-1);
    }

  lower_mass=atof(argv[1]); 
  ifstream ifs(argv[2]);

  double fs[number_z+1][100];
  double ms[100];
  double zs[number_z+1];
  double z,m,f1;
  double n=0;

  for(int i=0;i<=number_z;++i)
    {
  
   
      for(int j=0;j<100;++j)
	{
	  ifs>>z>>m>>f1;
	  zs[i]=z;
	  ms[j]=m;
	  fs[i][j]=f1;
	}
    
    }
  
  
  for(int i=0;i<number_z;++i)
    {
      n=0;
      for(int j=0;j<100;++j)
	{
	  if(ms[j]>=lower_mass&&zs[i]<=upper_redshift)
	    {
	      //double delta_z=zs[i+1]-zs[i];
	      double delta_m=ms[j]-ms[j-1];
	      //cout<<delta_z<<"\t"<<delta_m<<endl;
	      da=calc_angular_distance(zs[i])/Mpc;
	      n+=delta_m*fs[i][j]*(c/km)/(calc_H(zs[i])/(km/s/Mpc))*da*da*(1+zs[i])*(1+zs[i]);//*pow(10./180*3.1415926,2);
	    }
	}
      cout<<zs[i]<<"\t";
        cout<<n<<endl;
	
  
    }
  // cout<<lower_mass<<" "<<(int)n<<endl;

  //cout<<lower_mass<<" "<<n<<endl; 
}

