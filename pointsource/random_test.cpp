/*-
 * Test the distribution of random numbers generated
 * by TRNG library.
 *
 * Aaron LI
 * 2015/05/13
 */

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>

#include <trng/yarn3.hpp>
#include <trng/uniform_dist.hpp>
#include <trng/uniform01_dist.hpp>


int
main( int argc, char **argv )
{
    if ( argc != 3 ) {
        std::cerr << "Usage: " << std::endl
                  << "    " << argv[0] << " <number> <outfile>" << std::endl;
        exit( EXIT_FAILURE );
    }

    // number of points
    long number = atol( argv[1] );
    // output file
    std::ofstream outfile;
    outfile.open( argv[2] );

    // random number engine to be tested
    trng::yarn3 R;
    // initialize
    R.seed( time(NULL) );

    // random number distribution
    trng::uniform01_dist<> u01;
    trng::uniform_dist<> u55( -5.0, 5.0 );
    trng::uniform_dist<> u03( 0.0, 3.0 );

    double x1, x2, x3, y1, y2, y3;
    long i;

    outfile << "# x1\tx2\tx3\ty1\ty2\ty3" << std::endl;

    for ( i=0; i<number; i++ ) {
        x1 = u01( R );
        x2 = u55( R );
        x3 = u03( R );
        y1 = u01( R );
        y2 = u55( R );
        y3 = u03( R );
        outfile << x1 << "\t" << x2 << "\t" << x3 << "\t"
                << y1 << "\t" << y2 << "\t" << y3 << std::endl;
    }

    outfile.close();
    return 0;
}

