#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "adapt_trapezoid.h"
#include "calc_distance.h"

static const double cm=.01;
static const double m=100*cm;
static const double s=1;
static const double km=1000*100*cm;
static const double Mpc=3.08568e+24*cm;
static const double kpc=3.08568e+21*cm;
static const double yr=365.*24.*3600.;
static const double Gyr=1e9*yr;
static const double kg=1;
static const double N=kg*m/s/s;
static const double M_sun=1.99E30*kg;

static const double pi=3.141592653;
static const double omega_m=.27;
static const double omega_l=.73;
static const double H0=71*km/s/Mpc;//Hubble constant changed 70->71, 2009/08/20
static const double G=6.6726E-11*N*m*m/kg/kg;
const double h71=1;

static const double delta_c=200;

double ra,Dec,z,flux,rcDA,T,alpha,beta,mass;
using namespace std;

double foo(double r)
{
  //  double beta=2./3.;
  return 2*pi*r*pow(1+(r/rcDA)*(r/rcDA),-3.*beta+0.5);
}


double calc_H(double zz)
{
  return H0*E(zz);
}

double calc_rhoc(double zz)
{
  double H=calc_H(zz);
  return 3*H*H/8/pi/G;
}


int main(int argc,char* argv[])
{
  if (argc != 3)
    {
      cerr << "Usage:" << endl;
      cerr << "    " << argv[0] << " <in: ra_dec_z_m_F${nu}_rcDA_T_alpha.dat> "
           << "<out: ra_dec_F${nu}_S0_rcut_rcDA_beta.dat>" << endl;
      exit( EXIT_FAILURE );
    }

  ifstream ifs(argv[1]);
  ofstream ofs(argv[2]);

  double rho_c, r200, rcut, S0;

  while (1)
    { 
      ifs>>ra;
      if(!ifs.good())
	{
	  break;
	}
      ifs>>Dec>>z>>mass>>flux>>rcDA>>T>>alpha;
      
      rho_c=calc_rhoc(z);
      r200=pow(mass*M_sun/(rho_c*delta_c)*3./4./pi,1/3.);
      rcut=r200/(calc_angular_distance(z)/100.)/pi*180.*60.;
      beta=8.85E-15*mass/T/(r200*m/Mpc)/pow(rcut/rcDA,2)*(1+pow(rcut/rcDA,2));
      S0=flux/adapt_trapezoid(foo,0.,rcut,.001);
      //cerr<<adapt_trapezoid(foo,0.,rcut,.001)<<endl;
      //cerr<<foo(0.1)<<endl;;
      
      ofs << ra << " " << Dec << " " << flux << " " << S0 << " "
          << rcut << " " << rcDA << " " << beta <<endl;
    }

  return 0;
}  

