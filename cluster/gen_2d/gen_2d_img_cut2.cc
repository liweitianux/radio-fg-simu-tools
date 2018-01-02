/**
 * Make 2D emission image for all simulated radio halos,
 * as well as create a mask image to track the pixels of each halo.
 *
 * Based on 'gen_2d_img_cut.cc'
 *
 * WARNING:
 * Multiple halos may overlap with each other, how to deal with this??
 *
 * Weitian LI
 * 2017-10-16
 */

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <utility>
#include <fstream>
#include "fio.h"

#include <cassert>
#include <vector>

using namespace std;
using namespace ::blitz;
const double min_theta=85;

const double pi=3.14159265358979323846;
const double scale=10;
int width;

const int img_width=1024;
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
      cerr<<"Usage:"<<argv[0]<<" <data: ra_dec_Fv_s0_rcut_rcDA_beta.dat> <out: image> <out: mask>"<<endl;
      return -1;
    }


  std::vector<double> i_vec,j_vec,s0_vec,rc_vec;
  std::vector<double> rcut_vec,beta_vec;
  ifstream data_file(argv[1]);
  while(1)
    {
      double ra,dec,F_50_200,s0,rcut,rc,beta;
      data_file>>ra>>dec>>F_50_200>>s0>>rcut>>rc>>beta;
      if(!data_file.good())
	{
	  break;
	}
      pair<int,int> ij=sky2img(ra,dec);
      i_vec.push_back(ij.first);
      j_vec.push_back(ij.second);
      
      s0_vec.push_back(s0);

      
      rc_vec.push_back(rc/60./sky_width*img_width);
      rcut_vec.push_back(rcut/60./sky_width*img_width);
      beta_vec.push_back(beta);
      cerr<<ra<<"\t"<<dec<<"\t"<<ij.first<<"\t"<<ij.second<<"\t"<<*(rc_vec.end()-1)<<"\t"<<s0_vec.back()<<endl;
    }

  Array<double,2> mx(img_width,img_width);
  Array<double,2> mask(img_width,img_width);
  mx=0;
  mask=0;
  std::vector<double>::size_type k;
  for(k=0;k<i_vec.size();++k)
    {
      for(int i=-rcut_vec[k];i<rcut_vec[k];++i)
	{
	  for(int j=-rcut_vec[k];j<rcut_vec[k];++j)
	    {
	      double r=sqrt(i*i+j*j*1.);
	      if(i_vec[k]+i>=0&&
		 i_vec[k]+i<img_width&&
		 j_vec[k]+j>=0&&
		 j_vec[k]+j<img_width&&
		 r<rcut_vec[k])
		{
		  mx((int)(i+i_vec[k]),(int)(j+j_vec[k]))+=
		    beta_s(r,rc_vec[k],beta_vec[k],s0_vec[k]);
		  mask((int)(i+i_vec[k]),(int)(j+j_vec[k]))=k+1;
		}
	    }
	}
    }
  
  cfitsfile ff;
  ff.create(argv[2]);
  ff<<mx;

  cfitsfile fmask;
  fmask.create(argv[3]);
  fmask<<mask;
}
