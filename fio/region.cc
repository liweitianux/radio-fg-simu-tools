#include "region.h"

using namespace std;

region::region()
  :region_imp(0)
{}


region::~region()
{
  release();
}

region::region(const region& rhs)
{
  release();
  region_imp=(rhs.region_imp)->copy();
}

region::region(string str)
{
  region_imp=region_parser::parse(str);
  assert(region_imp!=0);
}

const region& region::operator=(const region& rhs)
{
  if(this!=&rhs)
    {
      release();
      region_imp=(rhs.region_imp)->copy();
    }
  return (*this);
}

void region::release()
{
  if(region_imp)
    {
      delete region_imp;
    }
}


string region::to_string()const
{
  assert(region_imp!=0);
  if(region_imp)
    {
      return region_imp->to_string();
    }
  else
    {
      return string("null");
    }
}

std::set<region::ipixel> region::get_pixel_set()const
{
  if(region_imp)
    {
      return region_imp->get_pixel_set();
    }
  else
    {
      return set<region::ipixel>();
    }
}

///end

