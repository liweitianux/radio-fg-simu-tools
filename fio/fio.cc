#include "fio.h"

fitsfile* openimage(const char* name)
{
  int status=0;
  fitsfile* ff;
  if(ffiopn(&ff,name,READWRITE,&status))
    {
      std::cerr<<"error:"<<status<<std::endl;
    }
  return ff;
}

fitsfile* createfitsfile(const char* name)
{
  int status=0;
  fitsfile* ff;
  remove(name);
  if(ffinit(&ff,name,&status))
    {
      std::cerr<<"error:"<<status<<std::endl;
    }
  return ff;
}

void closefitsfile(fitsfile*& ff)
{
  int status=0;
  if(fits_close_file(ff,&status))
    {
      std::cerr<<"error:"<<status<<std::endl;
      exit(-1);
    }
}






long get_num_rows(cfitsfile& ff)
{
  int status=0;
  long nrows=0;
  int hdutype;
  if(ffmahd(ff.fptr(),2,&hdutype,&status))
    {
      std::cerr<<"error:"<<status<<std::endl;
    }
  if(ffgnrw(ff.fptr(),&nrows,&status))
    {
      std::cerr<<"error!"<<status<<std::endl;
    }
  return nrows;
}




int get_col_num(cfitsfile& ff,const char* colname)
{
  int status=0;
  int colnum=0;
  if(ffgcno(ff.fptr(),CASEINSEN,const_cast<char*>(colname),&colnum,&status))
    {
      std::cerr<<"error:"<<status<<std::endl;
    }
  return colnum;
}


double get_x_offset(cfitsfile& ff)
{
  int status=0;
  double x_offset;
  if(fits_read_key(ff.fptr(),TDOUBLE,"LTV1",&x_offset,NULL,&status))
    {
      std::cerr<<status<<std::endl;
    }
  return x_offset;
}

double get_y_offset(cfitsfile& ff)
{
  int status=0;
  double y_offset;
  if(fits_read_key(ff.fptr(),TDOUBLE,"LTV2",&y_offset,NULL,&status))
    {
      std::cerr<<status<<std::endl;
    }
  return y_offset;
}

void put_offsets(cfitsfile& ff,double x_offset,double y_offset)
{
  int status=0;
  if(fits_update_key(ff.fptr(),TDOUBLE,"LTV1",&x_offset,NULL,&status))
    {
      std::cerr<<"error:"<<status<<std::endl;
    }
  if(fits_update_key(ff.fptr(),TDOUBLE,"LTV2",&y_offset,NULL,&status))
    {
      std::cerr<<"error:"<<status<<std::endl;
    }
  x_offset=-x_offset;
  y_offset=-y_offset;
  if(fits_update_key(ff.fptr(),TDOUBLE,"CRVAL1P",&x_offset,NULL,&status))
    {
      std::cerr<<"error:"<<status<<std::endl;
    }
  if(fits_update_key(ff.fptr(),TDOUBLE,"CRVAL2P",&y_offset,NULL,&status))
    {
      std::cerr<<"error:"<<status<<std::endl;
    }
}

void pass_offsets(cfitsfile& ff1,cfitsfile& ff2)
{
  double x_offset=get_x_offset(ff1);
  double y_offset=get_y_offset(ff1);
  put_offsets(ff2,x_offset,y_offset);
}

