#include "fio.h"

#define DOUBLE_MAX 9e300;
#define DOUBLE_MIN 0;
#define FLOAT_MAX 9e30;
#define FLOAT_MIN 0;


double fits_trait<double>::nulval=0;
float fits_trait<float>::nulval=0;
long fits_trait<long>::nulval=0;
short fits_trait<short>::nulval=0;

double fits_trait<double>::max=DOUBLE_MAX;
float fits_trait<float>::max=FLOAT_MAX;
long fits_trait<long>::max=LONG_MAX;
short fits_trait<short>::max=SHRT_MAX;

double fits_trait<double>::min=DOUBLE_MIN;
float fits_trait<float>::min=FLOAT_MIN;
long fits_trait<long>::min=LONG_MIN;
short fits_trait<short>::min=SHRT_MIN;
