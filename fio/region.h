#ifndef REGION_H
#define REGION_H
#include "region_imp.h"



class region
{
 private:
  region_base* region_imp;
 public:
  typedef region_base::pixel pixel;
  typedef region_base::ipixel ipixel;

 private:
  void release();
 public:
  region();
  ~region();
  region(const region&);
  region(std::string);


  const region& operator=(const region&);
  //bool operator==(const region&);
  //bool operator!=(const region&);

  
  std::string to_string()const;
  std::set<ipixel> get_pixel_set()const;
  
};


#endif

