#ifndef REGION_IMP_H
#define REGION_IMP_H
#include <string>
#include <set>
#include <sstream>
#include <cassert>
#include <iostream>

class region_base
{
 public:
  typedef std::pair<double,double> pixel;
  typedef std::pair<int,int> ipixel;
  virtual std::string to_string()const=0;
  virtual std::set<ipixel> get_pixel_set()const=0;
  virtual region_base* copy()const=0;

};


class circle_region
:public region_base
{
 public:
  typedef region_base::pixel pixel;
  typedef region_base::ipixel ipixel;
 public:
  pixel center;
  double radius;

  circle_region();
  circle_region(pixel c,double r);

  std::string to_string()const;
  std::set<ipixel> get_pixel_set()const;
  region_base* copy()const;
};


class pie_region
:public region_base
{
 public:
  typedef region_base::pixel pixel;
  typedef region_base::ipixel ipixel;
 public:
  pixel center;
  double radius1,radius2;
  double angule1,angule2;

  pie_region();
  pie_region(pixel c,double r1,double r2,double a1,double a2);

  std::string to_string()const;
  std::set<ipixel> get_pixel_set()const;
  region_base* copy()const;
};

class annulus_region
:public region_base
{
 public:
  typedef region_base::pixel pixel;
  typedef region_base::ipixel ipixel;
 public:
  pixel center;
  double radius1,radius2;
  
  annulus_region();
  annulus_region(pixel c,double r1,double r2);

  std::string to_string()const;
  std::set<ipixel> get_pixel_set()const;
  region_base* copy()const;
};


class region_parser
{
 public:
  typedef region_base::pixel pixel;
  static region_base* parse(std::string str)
    {
      if(str.find("circle")==0)
	{
	  str.erase(0,6);
	  str.erase(str.find("("),1);
	  str.replace(str.find(","),1," ");
	  str.replace(str.find(","),1," ");
	  
	  
	  std::istringstream iss(str);
	  double ii,jj,r;
	  iss>>jj>>ii>>r;
	  
	  region_base* ptr=new circle_region(pixel(ii,jj),r);
	  return ptr;
	}
      else if(str.find("pie")==0)
	{
	  str.erase(0,3);
	  str.erase(str.find("("),1);
	  str.replace(str.find(","),1," ");
	  str.replace(str.find(","),1," ");
	  str.replace(str.find(","),1," ");
	  str.replace(str.find(","),1," ");
	  str.replace(str.find(","),1," ");
	 
	  double ii,jj,r1,r2,a1,a2;
	  std::istringstream iss(str);
	  iss>>jj>>ii>>r1>>r2>>a1>>a2;
	  region_base* ptr=new pie_region(pixel(ii,jj),r1,r2,a1,a2);
	  return ptr;
	}
      else if(str.find("annulus")==0)
	{
	  str.erase(0,7);
	  str.erase(str.find("("),1);
	  str.replace(str.find(","),1," ");
	  str.replace(str.find(","),1," ");
	  str.replace(str.find(","),1," ");
	  double ii,jj,r1,r2;
	  std::istringstream iss(str);
	  iss>>jj>>ii>>r1>>r2;
	  region_base* ptr=new annulus_region(pixel(ii,jj),r1,r2);
	  return ptr;
	}
      return 0;
    }
};




#endif
