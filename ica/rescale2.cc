#include <blitz/array.h>
#include <fio.h>
#include <iostream>
#include <fstream>
#include <core/optimizer.hpp>
#include <methods/powell/powell_method.hpp>
#include <methods/gsl_simplex/gsl_simplex.hpp>
#include <vector>
#include <algorithm>
#include <sstream>

using namespace std;
using namespace blitz;
using namespace opt_utilities;


vector<blitz::Array<double,1> > imgs;
blitz::Array<double,1> I_target;
blitz::Array<short,2> mask;


blitz::Array<double,1> two_to_one(const blitz::Array<double,2>& rhs,
				  const blitz::Array<short,2>& mask)
{
  int size1=sum(mask);
  blitz::Array<double,1> result(size1);
  result=0;
  int n=0;
  for (int i=0;i<mask.extent(0);++i)
    {
      for(int j=0;j<mask.extent(1);++j)
	{
	  if(mask(i,j)!=0)
	    {
	      result(n)=rhs(i,j);
	      n+=1;
	    }
	}
    }
  return result;
}


double calc_bkg(blitz::Array<double,2>& img,
		blitz::Array<short,2>& mask,
		int center_i,int center_j,int r)
{
  double result=0;
  double cnt=0;
  double max_value=-1e99;;
  double min_value=1E99;
  for(double r1=r-1;r1<r+2;++r1)
    {
      for(double theta=0;theta<2.*3.1415926;theta+=.5/r1)
	{
	  double di=r1*cos(theta);
	  double dj=r1*sin(theta);
	  if(mask(int(di+center_i),int(dj+center_j))!=0)
	    {
	      result+=img(int(di+center_i),int(dj+center_j));
	      max_value=std::max(img(int(di+center_i),int(dj+center_j)),max_value);
	      min_value=std::min(img(int(di+center_i),int(dj+center_j)),min_value);
	      ++cnt;
	    }
	}
    }
  return min_value;
  return result/cnt;
}


blitz::Array<double,2> one_to_two(const blitz::Array<double,1>& rhs,
				  const blitz::Array<short,2>& mask)
{
  int size1=sum(mask);
  blitz::Array<double,2> result(mask.shape());
  result=0;
  int n=0;
  for(int i=0;i<mask.extent(0);++i)
    {
      for(int j=0;j<mask.extent(1);++j)
	{
	  if(mask(i,j)!=0)
	    {
	      result(i,j)=rhs(n);
	      n+=1;
	    }
	}
    }
  
  return result;
}


class foo
  :public func_obj<double,vector<double> >
{
public:
  double do_eval(const vector<double>& p)
  {
    blitz::Array<double,1> result(I_target.copy());
    result-=p[0];
    //vector<double> p1=p;
    //    p1[4]=10;
    for(int i=1;i<p.size();++i)
      {
	result-=p[i]*imgs[i-1];
      }
    
    //    result/=(I_target+0.01);
    double result1=sum(result*result);
    for(int i=0;i<p.size();++i)
      {
	cout<<p[i]<<" ";
      }
    cout<<result1<<"\n";
    return result1;

  }
  
  foo* do_clone()const
  {
    return new foo(*this);
  }
  
  
};


int main(int argc,char* argv[])
{

  if(argc<6)
    {

      cerr<<"Usage:"<<argv[0]<<" <mask> <target> <source1> [source2] ..."\
	"<source cluster> [output name]"<<endl;
      return -1;
    }

  
  cfitsfile ff_mask;
  ff_mask.open(argv[1]);
  ff_mask>>mask;

  int max_i=0;
  int min_i=65536;
  int max_j=0;
  int min_j=65535;

  for(int i=0;i<mask.extent(0);++i)
    {
      for(int j=0;j<mask.extent(1);++j)
	{
	  if(mask(i,j)>0)
	    {
	      max_i=std::max(i,max_i);
	      max_j=std::max(j,max_j);
	      min_i=std::min(i,min_i);
	      min_j=std::min(j,min_j);
	    }
	}
    }
  cout<<max_i<<endl;
  int center_i=(max_i+min_i)/2;
  int center_j=(max_j+min_j)/2;

  cout<<center_i<<"\t"<<center_j<<endl;
  
  int r=((max_i-min_i)+(max_j-min_j))/4;
  cout<<r<<endl;
  cfitsfile ff_target;


  Array<double,2> mx;
  
  ff_target.open(argv[2]);
  ff_target>>mx;
  
  Array<double,1> mx1(two_to_one(mx,mask));
  I_target.resize(mx1.shape());
  I_target=mx1;
  
  int freq=atoi(argv[1]);
  cout<<"sources:"<<endl;
  for(int i=3;i<=argc-2;++i)
    {
      
      cfitsfile ff;
      ff.open(argv[i]);
      cout<<argv[i]<<endl;
      //blitz::Array<double,2> img;
      ff>>mx;
      //imgs[i-1].resize(mx1.shape());
      imgs.push_back(two_to_one(mx,mask));
      
      //imgs[i-1]-=min(imgs[i-1]);
      double bkg=calc_bkg(mx,mask,center_i,center_j,r);
      imgs.back()-=bkg;
      cout<<bkg<<endl;
    }
  //  return 0;


  //  cout<<min(I_clu)<<endl;

  //  return 0;
 
  


  optimizer<double,vector<double> > opt;
  //t.set_func_obj(foo());
  foo f;
  opt.set_func_obj(f);
  //powell_method<double,vector<double> >pm;
  gsl_simplex<double,vector<double> > pm;

  opt.set_opt_method(pm);
  opt.set_precision(1E-10);

  vector<double> p(argc-4);

  opt.set_start_point(p);

  p=opt.optimize();

  cfitsfile ff_out;

  ff_out.create(argv[argc-1]);
  imgs.back()*=p.back();
  cout<< "fdfa";
  Array<double,2> outimg(one_to_two(imgs.back(),mask));
  ff_out<<outimg;
}
