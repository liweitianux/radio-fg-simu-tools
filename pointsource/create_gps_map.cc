#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <fio.h>

using namespace std;

const double pi = 3.141592654;
const double c = 2.99792458E8; // m/s
const double kb = 1.38E-23;
const double sky_size_degree = 10;
const int img_size = 1024;


double calc_Tb(double flux_in_Jy, double freq_in_MHz)
{
  double freq=freq_in_MHz*1E6;
  double Omegab=pow(sky_size_degree/img_size/180.0*pi, 2);
  double Sb=flux_in_Jy*1E-26/Omegab;
  double result=Sb/2.0/freq/freq*c*c/kb;
  return result;
}


int main(int argc, char **argv)
{
  if (argc != 4)
    {
      cerr << "Usage:" << endl;
      cerr << "    " << argv[0] << " <gps_data> <freq_MHz> <output>" << endl;
      exit( EXIT_FAILURE );
    }

  ifstream gps_data(argv[1]);
  double freq=atof(argv[2]);
  cfitsfile ff;
  ff.create(argv[3]);

  blitz::Array<double,2> mx(img_size,img_size);
  mx=0;
  double flux_5000,z,x,y,Dtrue,da_Mpc,S_vp,vp_GHz,k,l,a1,a2;
  int i = 0;

  for(;;)
    {
      gps_data>>flux_5000>>z>>x>>y
	      >>Dtrue>>da_Mpc
	      >>S_vp>>vp_GHz>>k>>l>>a1>>a2;
      if(!gps_data.good())
	{
	  break;
	}
      double v_GHz=freq/1000.;
      double Sv=S_vp*(1/(1-exp(-1))*pow(v_GHz/vp_GHz,k)*(1-exp(-pow(v_GHz/vp_GHz,l-k))));
      double Tb = calc_Tb(Sv,freq);
      int xi = (int) ((x+5)/sky_size_degree*img_size);
      int yi = (int) ((y+5)/sky_size_degree*img_size);
#ifdef DEBUG
      cerr << "DEBUG: " << i << "\t" << xi << "\t" << yi << "\t" << Tb << endl;
#endif /* DEBUG */
      if (xi >= img_size || yi >= img_size) {
          cerr << "WARNING: exceeded image size: xi=" << xi
               << ", yi=" << yi << endl;
      } else {
          mx(xi,yi) = Tb;
      }
      i++;
    }

  ff<<mx;

  return 0;
}

