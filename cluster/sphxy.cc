#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <time.h>
#include <cassert>

using namespace std;
const double pi=3.141592653;

int main(int argc,char* argv[])
{
  if (argc != 3)
    {
      cerr << "Usage:" << endl;
      cerr << "    " << argv[0] <<" <number> <out: ra_dec.dat>" << endl;
      exit(-1);
    }
  double n_max = atoi(argv[1]);
  ofstream ofs(argv[2]);
 
  srand(time(NULL));
 
  int n=1;
  double x, y, thi, phi;
 
  while(n<=n_max)
    {
      x = (double)rand()/RAND_MAX*10.-5.;
      y = (double)rand()/RAND_MAX*10.-5.;
    
      thi = 90.-sqrt(x*x+y*y);
      phi = atan2(x,y)/pi*180.+180.; 
      ofs << phi << "  " << thi << endl;
	  
      n++;
    }
}

