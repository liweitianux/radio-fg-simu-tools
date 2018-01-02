#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "f_rosat.h"

using namespace std;

double calc_Lx_bol(double Tx, double z)
{
  double Lx_bol = 12.44e44*pow(Tx/6.,2.64)*pow(1+z,1.52);
  return Lx_bol;
}

int main(int argc,char* argv[])
{
  if(argc != 3)
    {
      cerr << "Usage:" << endl;
      cerr << "    " << argv[0] <<" <in: z_m_T_P1400.dat> <out: z_m_T_P1400_Lrosat.dat>" << endl;
      exit(-1);
    }
  ifstream ifs(argv[1]);
  ofstream ofs(argv[2]);
  while(1)
    {
      double z,m,T,P1400;
      ifs>>z>>m>>T>>P1400;
      if(!ifs.good())
	{
	  break;
	}
      ofs << z << " " << m << " " << T << " " << P1400 << " "
          << calc_Lx_bol(T,z)/get_bolo_flux(z,T)*get_soft_flux(z,T) << endl;
    }

  return 0;
}

