//#define PGPLOT_USED

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <utility>
#include <fstream>
#include <map>
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

int img_width=2048;
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
  if(argc!=4)
    {
      cerr<<"Usage:"<<argv[0]<<" <data file ra dec s0 rc> <in image1> <in image2>"<<endl;
      return -1;
    }

#ifdef PGPLOT_USED
  pstream pout;
  pout.set_interactive(false);
#endif

  Array<double,2> mx1;
  Array<double,2> mx2;
  cfitsfile iff1;
  cfitsfile iff2;
  iff1.open(argv[2]);
  iff2.open(argv[3]);
  iff1>>mx1;
  iff2>>mx2;
  img_width=mx1.extent(0);
  
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
      

      std::map<pair<int,int>,double> flux_map1;
      std::map<pair<int,int>,double> flux_map2;
      for(double r=0;r<10*rc_vec.back();r+=.5)
	{
	  for(double theta=0;theta<3*pi;theta+=2*pi/(2*pi*r*5))
	    {
	      int ii=(int)(ij.first+r*cos(theta));
	      int jj=(int)(ij.second+r*sin(theta));

	      pair<int,int> p(ii,jj);
	      //	      cout<<p.first<<"\t"<<p.second<<"\t"<<mx1(ii,jj)<<endl;
	      
	      if(p.first>=0&&p.second>=0&&p.first<img_width&&p.second<img_width)
		{
		  double dflux1=mx1(ii,jj)*
		    ((60*sky_width)/img_width)*((60*sky_width)/img_width);
		  double dflux2=mx2(ii,jj)*
		    ((60*sky_width)/img_width)*((60*sky_width)/img_width);
		  //flux_map[p]=mx(ij.first,ij.second);//*((60*sky_width)/img_width)*((60*sky_width)/img_width);
		  flux_map1.insert(make_pair(p,
					     dflux1));
		  flux_map2.insert(make_pair(p,
					     dflux2));
		}
	    }
	}
      
      //cout<<flux_map.size()<<endl;
      double flux1=0;
      double flux2=0;
      for(map<pair<int,int>,double>::iterator i=flux_map1.begin();
	  i!=flux_map1.end();++i)
	{
	  //	  cout<<i->second<<endl;;
	  flux1+=i->second;
	}
      
      for(map<pair<int,int>,double>::iterator i=flux_map2.begin();
	  i!=flux_map2.end();++i)
	{
	  //	  cout<<i->second<<endl;;
	  flux2+=i->second;
	}

      cout<<ra<<"\t"<<dec<<"\t"<<flux1<<"\t"<<flux2<<"\t"<<flux1/flux2<<endl;
    }

  
  //  for(int i=img_width-1;i>=0;--i)
  
}
