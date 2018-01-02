#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

double data1,data2;
const double k=-2.1568;
const double b=3.45736;
double data3,data4;
double ra,Dec,flux408,s0,rcDA;

double beta_total=0;
int main(int argc,char* argv[])
{
  if (argc!=3)
    {
      cerr<<argv[0]<<" <z_m_T_P1400.dat> <ra_dec_F408_S0_rcDA.dat> "<<endl;
      return 0;
    }
  ifstream ifs(argv[1]);
  ifstream ifs2(argv[2]);
  while (1)  
    { 
      ifs>>data1;
      if (!ifs.good())
	{
	  break;
	}
      
      ifs>>data2>>data3>>data4;
      
      ifs2>>ra;
      if (!ifs2.good())
	{
	  break;
	}
      
      ifs2>>Dec>>flux408>>s0>>rcDA;
      double z=data1;
      double T=data3;

      double LL=12.44e44*pow((T/6.),2.64)*pow(1+z,1.52);
      // FIXME: 179 or 176/h71 ??
      double rrcc=179*pow(LL/5e44,0.2);
      double beta_out=0.00077*rrcc+0.472;
      cout<<beta_out<<endl;
      beta_total+=beta_out;
    }
  cout<<"!ave_beta="<<beta_total/1392<<endl;
}
