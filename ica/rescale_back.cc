#include <blitz/array.h>
#include <fio.h>
#include <iostream>
#include <fstream>
#include <core/optimizer.hpp>
#include <methods/powell/powell_method.hpp>
#include <methods/gsl_simplex/gsl_simplex.hpp>
#include <vector>
#include <algorithm>
#include <set>
#include <sstream>

using namespace std;
using namespace blitz;
using namespace opt_utilities;


vector<blitz::Array<double,1> > imgs;
blitz::Array<double,1> I_target;
blitz::Array<short,2> mask;
double max_bkg;
double sigma_bkg;

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
  double sum_value=0;
  double sum_value2=0;
  set<pair<int,int> > pix_set;

  for(double r1=r+1;r1<r+3;++r1)
    {
      for(double theta=0;theta<2.*3.1415926;theta+=.5/r1)
	{
	  double di=r1*cos(theta);
	  double dj=r1*sin(theta);
	  int ii=di+center_i;
	  int jj=dj+center_j;

	  set<pair<int,int> >::iterator iter=pix_set.find(make_pair(ii,jj));
	  
	  if(mask(int(di+center_i),int(dj+center_j))!=0&&iter==pix_set.end())
	    {
	      //	      cout<<"a"<<endl;
	      pix_set.insert(make_pair(ii,jj));
	      result+=img(int(di+center_i),int(dj+center_j));
	      max_value=std::max(img(int(di+center_i),int(dj+center_j)),max_value);
	      min_value=std::min(img(int(di+center_i),int(dj+center_j)),min_value);
	      ++cnt;
	      sum_value+=img(int(di+center_i),int(dj+center_j));
	      sum_value2+=img(int(di+center_i),int(dj+center_j))*
		img(int(di+center_i),int(dj+center_j));
	    }
	}
    }
  max_bkg=max_value;
  sigma_bkg=sqrt((sum_value2/cnt)-(sum_value/cnt)*(sum_value/cnt));
  //  return min(img);
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
    //    result-=imgs.back();
    
    //result/=(I_target);
    double result1=sum(result*result);
    for(int i=0;i<p.size();++i)
      {
	//	cout<<p[i]<<" ";
      }
    //    cout<<result1<<"\t"<<mean(result)<<"\n";
    return result1;

  }
  
  foo* do_clone()const
  {
    return new foo(*this);
  }
  
  
};


int is_negative(blitz::Array<double,2> img,int ci,int cj,int r)
{
  double sum1=0,sum2=0;
  int cnt1=0,cnt2=0;
  for(int i=0;i<img.extent(0);++i)
    {
      for(int j=0;j<img.extent(1);++j)
	{
	  double r1=sqrt((i-ci)*(i-ci)+(j-cj)*(j-cj));
	  if(r1<r/2)
	    {
	      sum1+=img(i,j);
	      ++cnt1;
	    }
	  else if(r1>=r/2 && r1<r)
	    {
	      sum2+=img(i,j);
	      ++cnt2;
	    }
	}
    }
  if (sum1/cnt1>=sum2/cnt2)
    {
      return false;
    }
  return true;
}


int main(int argc,char* argv[])
{

  if(argc<6)
    {

      cerr<<"Usage:"<<argv[0]<<" <mask> <input imgs> "\
	"<cluster> <target> <outfile>"<<endl;
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
  //cout<<max_i<<endl;
  int center_i=(max_i+min_i)/2;
  int center_j=(max_j+min_j)/2;

  //  cout<<center_i<<"\t"<<center_j<<endl;
  
  int r=((max_i-min_i)+(max_j-min_j))/4;
  
  cout<<"r="<<r<<endl;
  cfitsfile ff_target;

  string target_name(argv[argc-2]);

  Array<double,2> mx;
  
  ff_target.open(target_name.c_str());
  ff_target>>mx;
  
  Array<double,1> mx1(two_to_one(mx,mask));
  I_target.resize(mx1.shape());
  I_target=mx1;
  
  cout<<"------"<<endl;
  for(int i=2;i<=argc-3;++i)
    {
      cfitsfile ff;
      ff.open(argv[i]);
      //blitz::Array<double,2> img;
      ff>>mx;
      
      if(i==argc-3&&is_negative(mx,center_i,center_j,r))
	{
	  mx*=-1;
	  cout<<"negative\n";
	}



      //imgs[i-1].resize(mx1.shape());
      //imgs[i-1]=two_to_one(mx,mask);

      imgs.push_back(two_to_one(mx,mask));
      //imgs[i-1]-=min(imgs[i-1]);
      double bkg=calc_bkg(mx,mask,center_i,center_j,r);
      imgs.back()-=bkg;
      
      if(argc-3==i)
	{
	  cout<<argv[argc-3]<<endl;
	  cout<<bkg<<"\t"<<sigma_bkg<<"\t"<<
	    (sigma_bkg)/mean(imgs.back())<<endl;
	}
    }
  //  return 0;
  cout<<"------"<<endl;

  //  cout<<min(I_clu)<<endl;

  //  return 0;
 
  


  optimizer<double,vector<double> > opt;
  //t.set_func_obj(foo());
  foo f;
  opt.set_func_obj(f);
  powell_method<double,vector<double> >pm;
  //gsl_simplex<double,vector<double> > pm;

  opt.set_opt_method(pm);
  opt.set_precision(1E-10);

  vector<double> p(5);
  p[3]=1;
  opt.set_start_point(p);

  p=opt.optimize();

  cfitsfile ff_out;
  string outname(argv[argc-1]);
  

  cout<<p.back()<<endl;
  ff_out.create(outname.c_str());
  imgs.back()*=p.back();
  
  Array<double,2> outimg(one_to_two(imgs.back(),mask));
  ff_out<<outimg;
}
