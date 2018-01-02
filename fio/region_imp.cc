#include "region_imp.h"
#include <iostream>
#include <sstream>
#include <cmath>
using namespace std;

circle_region::circle_region()
{};

circle_region::circle_region(region_base::pixel c,double r)
  :center(c),radius(r)
{
  //  cout<<center.first<<endl;
  //  cout<<center.second<<endl;
  //  cout<<radius<<endl;
}

string circle_region::to_string()const
{

  ostringstream ostr;
  ostr<<"circle("<<center.second<<","<<center.first<<","<<radius<<")";
  return ostr.str();
}

set<circle_region::ipixel> circle_region::get_pixel_set()const
{
  const double pi=3.1415926;
  set<circle_region::ipixel> spixel;
  spixel.insert(circle_region::ipixel((int)(center.first),(int)(center.second)));
  for(double r=.5;r<=radius;r+=.3)
    {
      double angle_step=2*pi/(2*r*pi*2);
      for(double angle=0;angle<2*pi+.2;angle+=angle_step)
	{
	  int x,y;
	  x=static_cast<int>(center.second+r*cos(angle));
	  y=static_cast<int>(center.first+r*sin(angle));
	  spixel.insert(ipixel(y,x));
	}
    }
  return spixel;
}

region_base* circle_region::copy()const
{
  region_base* ptr =new circle_region(*this);
  return ptr;
}


pie_region::pie_region()
{};

pie_region::pie_region(region_base::pixel c,double r1,double r2,double a1,double a2)
  :center(c),radius1(r1),radius2(r2),angule1(a1),angule2(a2)
{
  //  cout<<center.first<<endl;
  //  cout<<center.second<<endl;
  //  cout<<radius<<endl;
}

string pie_region::to_string()const
{

  ostringstream ostr;
  ostr<<"pie("<<center.second<<","<<center.first<<","<<radius1<<","<<radius2<<","<<angule1<<","<<angule2<<")";
  return ostr.str();
}

set<pie_region::ipixel> pie_region::get_pixel_set()const
{
  const double pi=3.1415926;
  set<pie_region::ipixel> spixel;
  double r=radius1;
  if(r<1)
    {
      r=.5;
    }
  for(;r<=radius2;r+=.3)
    {
      double angle_step=2*pi/(2*r*pi*2);
      double a1=angule1/180*3.141592653589;
      double a2=angule2/180*3.141592653589;
      for(double angle=a1;angle<a2;angle+=angle_step)
	{
	  int x,y;
	  x=static_cast<int>(center.second+r*cos(angle));
	  y=static_cast<int>(center.first+r*sin(angle));
	  spixel.insert(ipixel(y,x));
	}
    }
  return spixel;
}

region_base* pie_region::copy()const
{
  region_base* ptr =new pie_region(*this);
  return ptr;
}
//////////////////////
annulus_region::annulus_region()
{};

annulus_region::annulus_region(region_base::pixel c,double r1,double r2)
  :center(c),radius1(r1),radius2(r2)
{
}

string annulus_region::to_string()const
{

  ostringstream ostr;
  ostr<<"annulus("<<center.second<<","<<center.first<<","<<radius1<<","<<radius2<<")";
  return ostr.str();
}

set<pie_region::ipixel> annulus_region::get_pixel_set()const
{
  const double pi=3.1415926;
  set<pie_region::ipixel> spixel;
  double r=radius1;
  if(r<1)
    {
      r=.5;
    }
  for(;r<=radius2;r+=.3)
    {
      double angle_step=2*pi/(2*r*pi*2);
      for(double angle=0;angle<2*3.14156+.1;angle+=angle_step)
	{
	  int x,y;
	  x=static_cast<int>(center.second+r*cos(angle));
	  y=static_cast<int>(center.first+r*sin(angle));
	  spixel.insert(ipixel(y,x));
	}
    }
  return spixel;
}

region_base* annulus_region::copy()const
{
  region_base* ptr =new annulus_region(*this);
  return ptr;
}


//end
