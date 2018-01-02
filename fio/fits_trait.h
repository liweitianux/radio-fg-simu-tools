#ifndef fits_traits_h
#define fits_traits_h
#include "fitsio.h"

template <typename T>
class fits_trait;
/*
{
 public:
  const static int bitpix=BYTE_IMG;
  const static int datatype=TBYTE;
};
*/
template <>
class fits_trait<short>
{
 public:
  const static int bitpix=SHORT_IMG;
  const static int datatype=TSHORT;
  static short nulval;
  static short max;
  static short min;
};


template <>
class fits_trait<long>
{
 public:
  const static int bitpix=LONG_IMG;
  const static int datatype=TLONG;
  static long nulval;
  static long max;
  static long min;
};

template <>
class fits_trait<float>
{
 public:
  const static int bitpix=FLOAT_IMG;
  const static int datatype=TFLOAT;
  static float nulval;
  static float max;
  static float min;
};

template <>
class fits_trait<double>
{
 public:
  const static int bitpix=DOUBLE_IMG;
  const static int datatype=TDOUBLE;
  static double nulval;
  static double max;
  static double min;
};


template <int typecode>
class typecode_trait;

template <>
class typecode_trait<TSHORT>
{
 public:
  typedef short datatype;
};

template <>
class typecode_trait<TLONG>
{
 public:
  typedef long datatype;
};

template <>
class typecode_trait<TFLOAT>
{
 public:
  typedef float datatype;
};

template <>
class typecode_trait<TDOUBLE>
{
 public:
  typedef double datatype;
};

/*
double fits_trait<double>::nulval=0;
float fits_trait<float>::nulval=0;
long fits_trait<long>::nulval=0;
short fits_trait<short>::nulval=0;
*/
#endif

