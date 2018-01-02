#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>

using namespace std;

const double h_70=1.;
double z,m;
double T;
double P_1400;
double omega_m=0.27;
double E(double z)
{
  return sqrt(omega_m*pow(1+z,3)+1-omega_m);
}

int main(int argc, char* argv[])
{ 
  if (argc != 3)
    {
      cerr << "Usage:" << endl;
      cerr << "    " << argv[0] << " <in: z_m.dat> <out: z_m_T_P1400.dat>" << endl;
      exit(-1);
    }

  ifstream ifs(argv[1]);
  ofstream ofs(argv[2]);

  while (1)
    {
      ifs>>z;
      if(!ifs.good())
	{
	  break;
	}
      ifs>>m;
      // T=(6.4*pow(0.70,2./3))*pow(m*1e-15,2./3)*pow(1+z,1);
      //T=pow(1.4*h_70*m/(5.75e13),1/1.6);//*pow(1+z,0.5/1.6);
      // T=pow(10.,(log((m/1e14)*E(z))/log(10)-0.74)/1.5)*6;       
      T=pow((E(z)*m/5.34e14),(1/1.72))*5;
      //T=pow((m/5.34e14),(1/1.72))*5;
      // P_1400=pow(10.,-0.390+9.83*(log(T/8)/log(10.)))*3.16e24/h_70; 
      //         P_1400=pow(10.,-0.24+6.40*(log(T/8)/log(10.)))*3.16e24/h_70;
      //      P_1400=pow(T,4.628)*2.8e20;
      P_1400=pow(T,2.88476)*1.072655e22*70*70/71./71.;

      ofs << z << " " << m << " " << T << " " << P_1400 << endl;
    }
}

