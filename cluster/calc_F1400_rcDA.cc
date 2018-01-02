#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <time.h>
#include <cassert>
#include "calc_distance.h"

using namespace std;

const double pi = 3.141592653;
const double h71 = 1;

int main(int argc, char* argv[])
{ 
  if (argc != 3)
    {
        cerr << "Usage:" << endl;
        cerr << "    " << argv[0] << "<in: z_m_T_P1400_Lrosat.dat> <out: z_m_F1400_rcDA_T.dat>" << endl;
        exit(-1);
    }
  ifstream ifs(argv[1]);
  ofstream ofs(argv[2]);

  double z,m,T,P_1400,Lrosat;
  // double alpha=atof(argv[1]);

  while (1)
    {
      ifs>>z;
      if(!ifs.good())
	{
	  break;
	}
   
      ifs>>m>>T>>P_1400>>Lrosat;
     
      double r_c=176/h71*pow((Lrosat/1e44/5),0.2);
      double d_a=calc_angular_distance(z);
      double rcDA=r_c*1000*3.086e18/d_a/pi*180.*60.;
      double d_l=calc_luminosity_distance(z);
      double flux_1400=P_1400*1e7/(4*pi*d_l*d_l);

      // double  k=flux_1400*pow(1400*1e6,alpha);
      //double F_408=pow(408./1400,-alpha)*flux_1400;
      //    double F_50_200=k/(1-alpha)*(pow(200*1e6,1-alpha)-pow(50*1e6,1-alpha));
      //ofs<<z<<" "<<T<<" "<<F_50_200<<" "<<rcDA<<endl;  
      ofs << z << " " << m << " " << flux_1400 << " " << rcDA
          << " " << T << endl;
    }
}

