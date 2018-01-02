//#define PGPLOT_USED

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <utility>
#include <fstream>
#include "fio.h"

#include <cassert>
#include <vector>
//#define PGPLOT_USED
#ifdef PGPLOT_USED
#include "pstream.h"
#include "ccpl_blitz.h"
using namespace ::ccpl;
using namespace ::ccpl::blitz;
#endif

using namespace std;
using namespace ::blitz;
const double min_theta=85;

const double pi=3.14159265358979323846;
const double scale=10;
int width;

const int img_width=2048;
const double sky_width=10;

int deg2pix(double deg)
{
  return (int)round(deg/sky_width*img_width);
}

double pix2deg(int pix)
{
  return round((double)pix/img_width*sky_width);
}


double beta_s(double r,double rc,double beta,double s0)
{
  double b=r/rc;
  return s0*pow(1+b*b,-3*beta+.5);
}

std::pair<int,int> sky2img(double ra,double dec)
{
  std::pair<int,int> result;
  result.first=(int)round(deg2pix(90.-dec)*(sin(ra/180.*pi))+img_width/2);
  result.second=(int)round(deg2pix(90.-dec)*(cos(ra/180.*pi))+img_width/2);
  return result;
}


int main(int argc,char* argv[])
{
  if(argc!=3)
    {
      cerr<<"Usage:"<<argv[0]<<" <data file ra dec flux s0 rcDA> <out image>"<<endl;
      return -1;
    }


#ifdef PGPLOT_USED
  pstream pout;
  pout.set_interactive(false);
#endif

  std::vector<double> i_vec,j_vec,s0_vec,rc_vec;
  ifstream data_file(argv[1]);
  while(1)
    {
      double ra,dec,F_50_200,s0,rc;
      data_file>>ra>>dec>>F_50_200>>s0>>rc;
      if(!data_file.good())
	{
	  break;
	}
      pair<int,int> ij=sky2img(ra,dec);
      i_vec.push_back(ij.first);
      j_vec.push_back(ij.second);
      
      s0_vec.push_back(s0);

      
      rc_vec.push_back(rc/60./sky_width*img_width);
      cout<<ra<<"\t"<<dec<<"\t"<<ij.first<<"\t"<<ij.second<<"\t"<<*(rc_vec.end()-1)<<endl;
    }

  Array<double,2> mx(img_width,img_width);
  //  for(int i=img_width-1;i>=0;--i)
  for(int i=0;i<img_width;++i)
    {
      for(int j=0;j<img_width;++j)
	{
	  double s=0;
          std::vector<double>::size_type k;
	  for(k=0;k<i_vec.size();++k)
	    {
	      double r=sqrt((i-i_vec[k])*(i-i_vec[k])+(j-j_vec[k])*(j-j_vec[k]));
	      //s+=exp(-d_degree/.1*d_degree/.1);
	      //s+=s0_vec[k]*pow(1+r*r/(rc_vec[k]*rc_vec[k]),-1.5);
	      double beta=0.66;
	      s+= beta_s(r,rc_vec[k],beta,s0_vec[k]);
	    }
	  mx(i,j)=s;
#if 0
	  if((i-img_width/2)*(i-img_width/2)+(j-img_width/2)*(j-img_width/2)<(img_width/2*img_width/2))
	    {
	      //	      mx(i,j)=s;
	    }
	  else
	    {
	      mx(i,j)=0;
	    }
#endif
	}
#ifdef PGPLOT_USED
      pout<<bitmap_plot(mx)<<endp;
#endif
      // cout<<i<<endl;
      cout<<(double)i/img_width*100<<"%"<<endl;
    }
  
  cout<<"fdsfsd"<<endl;
  cfitsfile ff;
  ff.create(argv[2]);
  ff<<mx;
}
