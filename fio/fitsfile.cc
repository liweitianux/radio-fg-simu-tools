#include "fitsfile.h"

cfitsfile::cfitsfile()
  :pff(NULL)
{}

cfitsfile::~cfitsfile()
{
  int status=0;
  if(pff)
    {
      if(fits_close_file(pff,&status))
	{
	  std::cerr<<"Error"<<status<<std::endl;
	}
    }
  pff=NULL;
}


fitsfile*& cfitsfile::fptr()
{
  return pff;
}

void cfitsfile::open(const char* name)
{
  int status=0;
  if(pff)
    {
      if(fits_close_file(pff,&status))
	{
	  std::cerr<<"Error:"<<status<<std::endl;
	}
    }
  if(ffopen(&pff,name,READWRITE,&status))
    {
      std::cerr<<"Error:"<<status<<std::endl;
    }
}


void cfitsfile::create(const char* name)
{
  int status=0;
  this->close();
  remove(name);  
  if(ffinit(&pff,name,&status))
    {
      std::cerr<<"Error:"<<status<<std::endl;
    }
}

void cfitsfile::close()
{
  int status=0;
  if(!pff)
    {
      return ;
    }
  if(fits_close_file(pff,&status))
    {
      std::cerr<<"Error"<<status<<std::endl;
    }
  pff=NULL;
}

