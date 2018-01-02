#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <set>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include "calc_distance.h"

using namespace std;

const double dark_mass=0.83;
double lower_mass;
int total_num;

double c=299792458;
double omega_E=0.73;//constant changed 2009/08/20
double omega_M=0.27;
double H_0=71;//Hubble constant changed 70->71 2009/08/20
const int number_z=152;
const int number_m=100;
static double Mpc=3.08567802e24;
static double km=1e5;
double n_z_m[number_z][number_m];
double n_z_lnm[number_z][number_m];
int Ncount_z_m[number_z][number_m];
double chi_square=0;

vector<double> z_line;
vector<double> lnm_line;

int main(int argc,char* argv[])
{
  if (argc != 5)
    {
      cerr << "Usage:" << endl;
      cerr << "    " << argv[0]
           << " DM_Mass_LowerLimit TotalNumber dndmdata <out: z_m.dat>" << endl;
      exit(-1);
    }

  lower_mass = atof(argv[1]);
  total_num = atoi(argv[2]);
  ifstream ifs(argv[3]);
  ofstream ofs(argv[4]);

  double z;
  double m;
  double dndm;
  srand(time(NULL));

  vector<double> N_z_m;
  set<double> z_set;
  set<double> m_set;
  vector<double> z_simu;
  vector<double> m_simu;
 
  while (1)
    {
      ifs>>z;
      if(!ifs.good())
	{
	  break;
	}
      ifs>>m;
      ifs>>dndm;

      z_set.insert(z);
      m_set.insert(m);
     
      double E_z=sqrt(omega_M*pow((1+z),3)+(1-omega_M-omega_E)*(1+z)*(1+z)+omega_E);
      double d_a=calc_angular_distance(z)/Mpc;
      double H_z=H_0*E_z;
      N_z_m.push_back(c/km/H_z*d_a*d_a*(1+z)*(1+z)*dndm*2*3.1415926*(1-cos(5./180*3.1415926)));

    }
  for (set<double>::iterator i=z_set.begin();i!=z_set.end();++i)
    {
      z_line.push_back(*i);
    }
  for (set<double>::iterator i=m_set.begin();i!=m_set.end();++i)
    {
      lnm_line.push_back(log(*i));
    }

  sort(z_line.begin(),z_line.end());
  sort(lnm_line.begin(),lnm_line.end());
  
  double max_n_z_lnm=0;
  for (unsigned int i=0;i<z_line.size();++i)
    {
      for (unsigned int j=0;j<lnm_line.size();++j)
	{
	  n_z_m[i][j]=N_z_m[number_m*i+j];

	  n_z_lnm[i][j]=exp(lnm_line[j])*n_z_m[i][j];
	  if (n_z_lnm[i][j]>max_n_z_lnm&&exp(lnm_line[j])>=lower_mass)
	    {
	      max_n_z_lnm=n_z_lnm[i][j];
	    }
	}
    }
  for (int i=1;i<number_z;++i)
    {
      for (int j=1;j<number_m+1;++j)
	{
	  Ncount_z_m[i][j]=0;
	}
    }
  ////////////////////////////////////////////////////////////
  //  double sum=0;
  //double count=0;
    
  for (int cnt=0;cnt<total_num;)
    {
      double z_rand=((double)rand())/RAND_MAX*(z_line.back()-z_line[0])+z_line[0];
      double lnm_rand=((double)rand())/RAND_MAX*(lnm_line.back()-log(lower_mass))+log(lower_mass);
      int z_upper=distance(z_line.begin(),upper_bound(z_line.begin(),z_line.end(),z_rand));
      int lnm_upper=distance(lnm_line.begin(),upper_bound(lnm_line.begin(),lnm_line.end(),lnm_rand));
      double t=(z_rand-z_line[z_upper-1])/(z_line[z_upper]-z_line[z_upper-1]);
      double u=(lnm_rand-lnm_line[lnm_upper-1])/(lnm_line[lnm_upper]-lnm_line[lnm_upper-1]);
      double y1=n_z_lnm[z_upper-1][lnm_upper-1];
      double y2=n_z_lnm[z_upper][lnm_upper-1];
      double y3=n_z_lnm[z_upper][lnm_upper];
      double y4=n_z_lnm[z_upper-1][lnm_upper];

      double cal_n_z_lnm=(1-t)*(1-u)*y1+t*(1-u)*y2+t*u*y3+(1-t)*u*y4;
      double rand_n=((double)rand())/RAND_MAX*max_n_z_lnm;

      //      count+=1;
      //sum+=cal_n_z_lnm*exp(lnm_rand);

      if (rand_n<cal_n_z_lnm)
	{
	  if (exp(lnm_rand)>lower_mass)
	    {
	     
	      double Nthero_z_m=n_z_m[z_upper][lnm_upper]*(z_line[z_upper]-z_line[z_upper-1])*(exp(lnm_line[lnm_upper])-exp(lnm_line[lnm_upper-1]));
#ifdef WANG_2010
//this expression must be included for simulating Wang et al. 2010 clusters.
	      if (Ncount_z_m[z_upper][lnm_upper]<=Nthero_z_m)
#endif
		{
		  ++Ncount_z_m[z_upper][lnm_upper];
		  chi_square+=pow((Ncount_z_m[z_upper][lnm_upper]-Nthero_z_m),2)/Nthero_z_m;
		  //	  ofs<<Ncount_z_m[z_upper][lnm_upper]<<"\t"<<Nthero_z_m<<"\t"<<chi_square<<endl;
		  ofs<<z_rand<<" "<<exp(lnm_rand)/dark_mass<<endl;
		  ++cnt;
		}
	      /*   else {
		   ofs<<Ncount_z_m[z_upper][lnm_upper]<<"\t"<<Nthero_z_m<<endl;

		   }*/
	    }
	}
    }
}

