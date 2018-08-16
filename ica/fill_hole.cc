#include <fio.h>
#include <cmath>
#include <iostream>

using namespace std;
using namespace blitz;


const int img_size=1024;

double search_nn_value(Array<double,2>& img,Array<short,2>& mask,
		       int ii,int jj)
{
  for(int r=1;r<20;++r)
    {
      double sum=0;
      double cnt=0;

      for(double theta=0;theta<2*3.1415926;theta+=.5/r)
	{
	  int di=r*cos(theta);
	  int dj=r*sin(theta);
	  
	  if(di+ii>=0&&dj+jj>=0&&di+ii<mask.extent(0)&&dj+jj<mask.extent(1))
	    {
	      if(mask(di+ii,dj+jj)>0)
		{
		  cnt+=1;
		  sum+=img(di+ii,dj+jj);
		}
	    }
	}

      if(cnt>4)
	{
	  return sum/cnt;
	}
    }
  return 0;
}





int main(int argc,char* argv[])
{
  if(argc!=4)
    {
      cerr<<"Usage:"<<argv[0]<<" <input image> <mask> <out image>"<<endl;
      return -1;
    }

  Array<double,2> in_img,out_img;
  Array<short,2> mask;
  cfitsfile ff;
  ff.open(argv[1]);
  ff>>in_img;

  cfitsfile ff_mask;
  ff_mask.open(argv[2]);
  ff_mask>>mask;

  out_img.resize(mask.shape());


  for(int i=0;i<mask.extent(0);++i)
    {
      cout<<i<<endl;
      for(int j=0;j<mask.extent(1);++j)
	{
	  if(mask(i,j)==0&&(i-img_size/2)*(i-img_size/2)+(j-img_size/2)*(j-img_size/2)<img_size*img_size/4*1.3)
	    //if(mask(i,j)==0)
	    {
	      double a=search_nn_value(in_img,mask,i,j);
	      //	      cout<<a<<endl;
	      out_img(i,j)=a;
	    }
	}
    }




  cfitsfile ff_out;
  ff_out.create(argv[3]);
  out_img+=(in_img*mask);
  ff_out<<out_img;
}
