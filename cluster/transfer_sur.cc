#include <iostream>
#include <cmath>
#include <cassert>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <fio.h>

using namespace std;
using namespace blitz;

const double c = 2.99792458e10; /* speed of light (cm/s) */
const double kb = 1.38e-16; /* Boltzmann constant (erg/k) */

const int img_width = 1024;

int
main(int argc, char **argv)
{
    if ( argc != 4 ) {
        cerr << "Usage:" << endl;
        cerr << "    " << argv[0] <<" <freq_MHz> <in: erg_img> <out: T_img>"
             << endl;
        exit( EXIT_FAILURE );
    }

    double nu = atof(argv[1]) * 1.0e6; /* MHz -> Hz */
    cfitsfile erg_img;
    erg_img.open(argv[2]);
    cfitsfile T_img;
    T_img.create(argv[3]);

    blitz::Array<double,2> erg_img_data;
    erg_img >> erg_img_data;

    blitz::Array<double,2> T_img_data(img_width, img_width);
    for (int i=0; i<img_width; ++i) {
        for (int j=0; j<img_width; ++j) {
            T_img_data(i, j) = erg_img_data(i, j) /
                (8.4681e-8 * 2.0 * kb * nu * nu / c / c);
        } 
    }
    T_img << T_img_data;

    return 0;
}
  
