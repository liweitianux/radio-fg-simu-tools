#include <iostream>
#include <cstdlib>
#include <fstream>
#include <random/normal.h>

using namespace std;
using namespace ranlib;

int
main(int argc, char* argv[])
{
    if (argc != 3) {
        cerr << "Usage:" << endl;
        cerr << "    " << argv[0] << " <number> <out: alpha.dat>" << endl;
        exit(-1);
    }

    int total_number = atoi(argv[1]);
    ofstream ofs(argv[2]);

    double alpha;
    double alpha_mean = 1.3;
    double alpha_std = 0.3;
    Normal<double>(0,1).seed((unsigned int)time(NULL));
    Normal<double> rnd(alpha_mean, alpha_std);
 
    int n = 0;
    while (n < total_number) {
        alpha = rnd.random();
        if (alpha > 1.0) {
            ofs << alpha << endl;
            n++;
        }
    }

    return 0;
}

