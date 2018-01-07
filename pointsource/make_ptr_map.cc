/*
 * Copyright (c) 2010 Junhua GU, Jingying Wang
 * Copyright (c) 2017-2018 Weitian LI
 *
 *
 * Simulate maps of point sources based on the database by Wilman et al. 2008
 * http://s-cubed.physics.ox.ac.uk/s3_sex/database_structure
 *
 * NOTE: The simulating sky size is fixed to be 10x10 [deg^2], which is
 *       constrained by the obtained database from Wilman et al. 2008.
 *       The obtained database covers R.A. -5:5 [deg] and Dec. -5:5 [deg].
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <cassert>
#include <string>
#include <vector>
#include <memory>
#include <cmath>
#include <CCfits/CCfits>

#include <unistd.h>  // getopt()

#include "matrix.h"


using namespace std;


const int FRI = 1;
const int FRII = 2;
const int RQQ = 3;
const int SF = 4;
const int SB = 5;
const int simulate_switches[] = {
        -1,  /* holder */
        1,   /* FR-I */
        1,   /* FR-II */
        1,   /* RQQ */
        1,   /* SF */
        1,   /* SB */
};


struct FRI_info {
    double core_ra;
    double core_dec;
    double core_flux;
    double lobe_flux[2];
    double lobe_major_axis[2];
    double lobe_minor_axis[2];
    double lobe_ra[2];
    double lobe_dec[2];
    double lobe_pa[2];
    double redshift;
    double cosv;

    void reset(void) {
        core_ra = -1;
        core_dec = -1;
        core_flux = -1;
        for (int i = 0; i < 2; ++i) {
            lobe_flux[i] = -1;
            lobe_major_axis[i] = -1;
            lobe_minor_axis[i] = -1;
            lobe_ra[i] = -1;
            lobe_dec[i] = -1;
            lobe_pa[i] = -1;
        }
        redshift = -1;
        cosv = -1;
    }
};

struct FRII_info {
    double core_ra;
    double core_dec;
    double core_flux;
    double lobe_flux[2];
    double lobe_major_axis[2];
    double lobe_minor_axis[2];
    double lobe_ra[2];
    double lobe_dec[2];
    double lobe_pa[2];
    double hotspot_flux[2];
    double hotspot_ra[2];
    double hotspot_dec[2];
    double cosv;
    double redshift;

    void reset(void) {
        core_ra = -1;
        core_dec = -1;
        core_flux = -1;
        for (int i = 0; i < 2; ++i) {
            lobe_flux[i] = -1;
            lobe_major_axis[i] = -1;
            lobe_minor_axis[i] = -1;
            lobe_ra[i] = -1;
            lobe_dec[i] = -1;
            lobe_pa[i] = -1;
            hotspot_flux[i] = -1;
            hotspot_ra[i] = -1;
            hotspot_dec[i] = -1;
            cosv = -1;
            redshift = -1;
        }
    }
};

struct RQQ_info {
    double flux;
    double ra;
    double dec;
    double redshift;

    void reset(void) {
        flux = -1;
        ra = -1;
        dec = -1;
        redshift = -1;
    }
};

struct SF_info {
    double flux;
    double ra;
    double dec;
    double major_axis;
    double minor_axis;
    double pa;
    double m_hi;
    int type;
    double redshift;

    void reset(void) {
        flux = -1;
        ra = -1;
        dec = -1;
        major_axis = -1;
        minor_axis = -1;
        pa = -1;
        m_hi = -1;
        redshift = -1;
    }
};


double rqq_spec(double I_151, double freqMHz)
{
    return pow(freqMHz/151, -0.7)*I_151;
}

double fr1_lobe_spec(double I_151, double freqMHz)
{
    return pow(freqMHz/151, -0.75)*I_151;
}

double fr1_core_spec(double I_151, double freqMHz)
{
    double freq = freqMHz * 1e6;  // [Hz]
    double a0 = (log10(I_151) - 0.07*log10(151E6) +
                 0.29*log10(151E6)*log10(151E6));
    double lgS = (a0 + 0.07*log10(freq) -
                  0.29*log10(freq)*log10(freq));
    return pow(10.0, lgS);
}

double fr2_core_spec(double I_151, double freqMHz)
{
    double freq = freqMHz * 1e6;  // [Hz]
    double a0 = (log10(I_151) - 0.07*log10(151E6) +
                 0.29*log10(151E6)*log10(151E6));
    double lgS = (a0 + 0.07*log10(freq) -
                  0.29*log10(freq)*log10(freq));
    return pow(10.0, lgS);
}

double fr2_lobe_spec(double I_151, double freqMHz)
{
    return pow(freqMHz/151, -0.75)*I_151;
}

double fr2_hotspot_spec(double I_151, double freqMHz)
{
    return pow(freqMHz/151, -0.75)*I_151;
}

double sf_spec(double I_151, double freqMHz)
{
    return pow(freqMHz/151, -0.7)*I_151;
}

double sb_spec(double I_151,double freqMHz)
{
    return pow(freqMHz/151, -0.7)*I_151;
}


double rad2deg(double rad)
{
    return rad * 180.0 / M_PI;
}

double deg2rad(double deg)
{
    return deg * M_PI / 180.0;
}


/*
 * flux -> brightness temperature [K]
 *
 * flux: [Jy]
 * omega: [arcsec^2]
 * freq: [MHz]
 */
double calc_Tb(double flux, double omega, double freq)
{
    const double c = 2.99792458e8;  // [m/s]
    const double k_B = 1.38e-23;
    const double arcsec2rad = M_PI / 180.0 / 3600.0;
    freq *= 1e6;  // [MHz]->[Hz]
    omega *= (arcsec2rad*arcsec2rad);
    double Sb = flux*1e-26 / omega;
    return 0.5 * Sb * c*c / (freq*freq * k_B);
}


/*
 * Determine whether the point (x, y) is inside the ellipse
 * defined with following properties:
 *    xc, yc : ellipse center
 *    a : semi-major axis
 *    b : semi-minor axis
 *    pa : position angle
 *    c = sqrt(a*a - b*b)
 *    f1x = xc + c*sin(pa)
 *    f1y = yc + c*cos(pa)
 *    f2x = xc - c*sin(pa)
 *    f2y = yc - c*cos(pa)
 *
 * XXX: why this calculation method?
 */
inline bool in_ellipse(int x, int y,
                       double f1x, double f1y,
                       double f2x, double f2y,
                       double a)
{
    double t1 = sqrt((x-f1x)*(x-f1x) + (y-f1y)*(y-f1y));
    double t2 = sqrt((x-f2x)*(x-f2x) + (y-f2y)*(y-f2y));
    return (t1+t2 < 2*a);
}


void usage(const char *name) {
    fprintf(stderr, "Usage: %s [ -o output_prefix ] "
            "[ -f fov ] [ -s image_size ] [ -x reduction_factor ] "
            "[ -m flux_min ] [ -M flux_max ] "
            "-O <output.csv> -i <dbfiles.list> freqMHz ...\n", name);
    fprintf(stderr, "\n");
    fprintf(stderr, "    -o : output prefix; default: 'ptr_'\n");
    fprintf(stderr, "    -f : sky size/fov [deg]; default: 10 [deg]\n");
    fprintf(stderr, "    -s : image size; default: 1800\n");
    fprintf(stderr, "    -x : reduce source number by the given factor\n");
    fprintf(stderr, "    -m : [uJy] minimum flux limit\n");
    fprintf(stderr, "    -M : [uJy] maximum flux limit\n");
    fprintf(stderr, "    -O : output CSV file storing the simulated source properties [required]\n");
    fprintf(stderr, "    -i : filename of list of input database files; [required]\n");
    fprintf(stderr, "    freqMHz : list of simulating frequencies [MHz]; [required]\n");
    fprintf(stderr, "\n");
    exit(EXIT_FAILURE);
}


int main(int argc, char * const argv[])
{
    /*
    * print command
    */
    cerr << "COMMAND:" << endl;
    for (auto && str : std::vector<std::string> { argv, argv + argc }) {
        cerr << str << " ";
    }
    cerr << endl << endl;

    const char *progname = argv[0];
    int opt;
    int img_size = 1800;
    int reduction_factor = 1;
    double fov_deg = 10.0;  // 10 [deg]
    double fov_asec = fov_deg * 3600.0;  // [arcsec]
    double flux_min = 0.0;  // [Jy]
    double flux_max = std::numeric_limits<double>::infinity();
    const char *output_prefix = "ptr_";
    const char *output_csv = NULL;
    const char *dbfiles_list = NULL;

    while ((opt = getopt(argc, argv, "f:i:m:o:s:x:M:O:")) != -1) {
        switch (opt) {
        case 'f':
            fov_deg = atof(optarg);  // [deg]
            fov_asec = fov_deg * 3600.0;  // [arcsec]
            break;
        case 'i':
            dbfiles_list = optarg;
            break;
        case 'm':
            flux_min = atof(optarg) / 1e6;  // [uJy]->[Jy]
            break;
        case 'o':
            output_prefix = optarg;
            break;
        case 's':
            img_size = atoi(optarg);
            break;
        case 'x':
            reduction_factor = atoi(optarg);
            break;
        case 'M':
            flux_max = atof(optarg) / 1e6;  // [uJy]->[Jy]
            break;
        case 'O':
            output_csv = optarg;
            break;
        default:  // '?'
            usage(progname);
        }
    }

    argc -= optind;
    argv += optind;
    if (argc == 0) {
        usage(progname);
    }
    if (output_csv == NULL) {
        usage(progname);
    }
    if (dbfiles_list == NULL) {
        usage(progname);
    }

    const double pix_deg = fov_deg / img_size;  // [deg]
    const double pix_size = fov_asec / img_size;  // [arcsec]
    const double pix_area = pix_size * pix_size;  // [arcsec^2]
    const double ra_min = -fov_deg / 2.0;
    const double dec_min = -fov_deg / 2.0;
    printf("INFO: flux limit: [%.2g, %.2g] [Jy]\n", flux_min, flux_max);
    printf("INFO: reduction factor: %d\n", reduction_factor);
    printf("INFO: output source data file: %s\n", output_csv);
    printf("INFO: input database files list: %s\n", dbfiles_list);
    printf("INFO: image size: %d\n", img_size);
    printf("INFO: pixel size: %.1f [arcsec]\n", pix_size);
    printf("INFO: sky fov: %.1f [deg]\n", fov_deg);

    vector<double> freqs;
    vector<Matrix<double> > images;
    printf("INFO: frequencies [MHz]:");
    for (int i = 0; i < argc; ++i) {
        freqs.push_back(atof(argv[i]));
        printf(" %.2f", freqs.back());
        images.push_back(Matrix<double>(img_size, img_size));
    }
    printf("\nINFO: number of frequencies: %lu\n", freqs.size());

    /*
    * Simulation
    */
    string fname;
    int old_galaxy;
    int lobe_cnt;
    int hotspot_cnt;
    int source_type;

    /* Database columns */
    long source,      /* 1. */
         cluster,     /* 2. */
         galaxy,      /* 3. */
         sftype,      /* 4. */
         agntype,     /* 5. */
         structure;   /* 6. */
    double ra,                /* 7. [deg] */
           dec,               /* 8. [deg] */
           distance,          /* 9. [cMpc]*/
           redshift,          /* 10. */
           pa,                /* 11. [rad] */
           major_axis,        /* 12. [arcsec] */
           minor_axis,        /* 13. [arcsec] */
           i_151,             /* 14. log10[Jy] */
           i_610,             /* 15. log10[Jy] */
           i_1400,            /* 16. log10[Jy] */
           i_4860,            /* 17. log10[Jy] */
           i_18000,           /* 18. log10[Jy] */
           m_hi,              /* 19. log10[Msun] */
           cos_va;            /* 20. cos(viewing_angle) */
    double flux;  // [Jy]

    ifstream input_list(dbfiles_list);
    ofstream csvout(output_csv);
    csvout << "# type: 1 - FRI, 2 - FRII, 3 - RQQ, 4 - SF, 5 - SB\n"
           << "# redshift\n"
           << "# x, y: (core) position on the simulated image\n"
           << "# flux: (core) flux @151[MHz]; unit: [Jy]\n"
           << "# major, minor, pa: (core) major and minor axes (unit: pixel), and position angle (unit: deg) (type: SF/SB)\n"
           << "# x_l1, y_l1, flux_l1, major_l1, minor_l1, pa_l1: lobe1 properties (type: FR-I/FR-II)\n"
           << "# x_l2, y_l2, flux_l2, major_l2, minor_l2, pa_l2: lobe2 properties (type: FR-I/FR-II)\n"
           << "# x_h1, y_h1, flux_h1: hotspot1 properties (type: FR-II)\n"
           << "# x_h2, y_h2, flux_h2: hotspot2 properties (type: FR-II)\n"
           << "#\n"
           << "type,redshift,x,y,flux,major,minor,pa"
           << ",x_l1,y_l1,flux_l1,major_l1,minor_l1,pa_l1"
           << ",x_l2,y_l2,flux_l2,major_l2,minor_l2,pa_l2"
           << ",x_h1,y_h1,flux_h2,x_h2,y_h2,flux_h2"
           << endl;

    long ptr_cnt_sum = 0;
    for (;;) {
        getline(input_list, fname);
        if (!input_list.good())
            break;

        cout << "Database: " << fname << endl;
        ifstream ifs(fname.c_str());

        old_galaxy = -1;
        lobe_cnt = 0;
        hotspot_cnt = 0;

        RQQ_info rqq;
        FRI_info fr1;
        FRII_info fr2;
        SF_info sf;

        long n = 0;
        long ptr_cnt = 0;
        for (n = 0; ; ++n) {
            ifs >> source;
            if (!ifs.good()) {
                cout << "  Source/component counts: " << n << endl;
                cout << "  Selected point sources: " << ptr_cnt << endl;
                break;
            }

            ifs >> cluster >> galaxy >> sftype >> agntype >> structure
                >> ra >> dec >> distance >> redshift
                >> pa >> major_axis >> minor_axis
                >> i_151 >> i_610 >> i_1400 >> i_4860 >> i_18000;
            if (sftype) {
                ifs >> m_hi;
            }
            ifs >> cos_va;

            flux = pow(10.0, i_151);  // [Jy]

            /*
            * ignore by galaxy, no by structure/component
            */
            if (old_galaxy != galaxy) {
                hotspot_cnt = 0;
                lobe_cnt = 0;
                fr1.reset();
                fr2.reset();
                rqq.reset();
                sf.reset();

                if (n % reduction_factor != 0)
                    continue;

                if (flux < flux_min || flux > flux_max)
                    continue;
            }

            if (sftype == 1)
                source_type = SF;
            else if (sftype == 2)
                source_type = SB;
            else if (agntype == 1)
                source_type = RQQ;
            else if (agntype == 2)
                source_type = FRI;
            else if (agntype == 3)
                source_type = FRII;
            else {
                cerr << "SF_type: " << sftype
                     << ", AGN_type: " << agntype << endl;
                assert(0);
            }

            /*
            * Radio-quiet AGN
            */
            if (source_type == RQQ) {
                if (simulate_switches[RQQ]) {
                    rqq.reset();
                    rqq.ra = ra - ra_min;
                    rqq.dec = dec - dec_min;
                    rqq.flux = flux;
                    rqq.redshift = redshift;

                    int x = (int)round(rqq.ra/pix_deg);
                    int y = (int)round(rqq.dec/pix_deg);
                    if (x < 0 || x >= img_size || y < 0 || y >= img_size)
                        continue;

                    for (vector<double>::size_type k = 0;
                         k < freqs.size();
                         ++k) {
                        images[k](x, y) += calc_Tb(
                            rqq_spec(rqq.flux, freqs[k]),
                            pix_area, freqs[k]);
                    }

                    ptr_cnt++;
                    csvout << source_type << "," << redshift << ","
                           << x << "," << y << "," << rqq.flux << endl;
                }
            }

            /*
            * Star-formation / starburst galaxies
            */
            if (source_type == SF || source_type == SB) {
                assert(structure == 4);
                sf.reset();
                sf.m_hi = m_hi;
                sf.ra = ra - ra_min;
                sf.dec = dec - dec_min;
                sf.flux = flux;
                sf.major_axis = major_axis;
                sf.minor_axis = minor_axis;
                sf.pa = pa;
                sf.redshift = redshift;
                sf.type = sftype;

                if ((simulate_switches[SF] && source_type == SF) ||
                    (simulate_switches[SB] && source_type == SB)) {
                    double x = sf.ra / pix_deg;
                    double y = sf.dec / pix_deg;
                    if (x < 0 || x >= img_size || y < 0 || y >= img_size)
                      continue;

                    double a = 0.5 * sf.major_axis / pix_size;
                    double b = 0.5 * sf.minor_axis / pix_size;
                    double c = sqrt(a*a - b*b);
                    double f1x = x + c*sin(pa);
                    double f1y = y + c*cos(pa);
                    double f2x = x - c*sin(pa);
                    double f2y = y - c*cos(pa);

                    double s = a*b * M_PI*pix_area;  // [arcsec^2]
                    if (s < pix_area)
                        s = pix_area;

                    int xmin = (int) round(x - a);
                    int xmax = (int) round(x + a);
                    int ymin = (int) round(y - a);
                    int ymax = (int) round(y + a);
                    for (int i = xmin; i <= xmax; ++i) {
                        for (int j = ymin; j <= ymax; ++j) {
                            if (i < 0 || i >= img_size || j < 0 || j >= img_size)
                                continue;

                            if (in_ellipse(i, j, f1x, f1y, f2x, f2y, a)) {
                                if (sf.type == 1) {
                                    /* star-forming */
                                    for (vector<double>::size_type k = 0;
                                         k < freqs.size();
                                         ++k) {
                                        images[k](i, j) += calc_Tb(
                                            sf_spec(sf.flux, freqs[k]),
                                            s, freqs[k]);
                                    }
                                } else if (sf.type==2) {
                                    /* starburst */
                                    for (vector<double>::size_type k = 0;
                                         k < freqs.size();
                                         ++k) {
                                        images[k](i, j) += calc_Tb(
                                            sb_spec(sf.flux, freqs[k]),
                                            s, freqs[k]);
                                    }
                                } else {
                                    assert(0);
                                }
                            }
                        }
                    }

                    ptr_cnt++;
                    csvout << source_type << "," << redshift << ","
                           << x << "," << y << "," << sf.flux << ","
                           << a*2 << "," << b*2 << ","
                           << rad2deg(pa) << endl;
                }
            }

          /*
           * FR II AGN
           */
            if (source_type == FRII) {
                assert(structure >= 1 && structure <= 3);
                if (structure == 1) {
                    // core
                    fr2.reset();
                    fr2.core_ra = ra-ra_min;
                    fr2.core_dec = dec-dec_min;
                    fr2.core_flux = flux;
                    fr2.redshift = redshift;
                    fr2.cosv = cos_va;
                } else if (structure == 2) {
                    // lobe
                    assert(lobe_cnt < 2);
                    fr2.lobe_flux[lobe_cnt] = flux;
                    fr2.lobe_ra[lobe_cnt] = ra - ra_min;
                    fr2.lobe_dec[lobe_cnt] = dec - dec_min;
                    fr2.lobe_major_axis[lobe_cnt] = major_axis;
                    fr2.lobe_minor_axis[lobe_cnt] = minor_axis;
                    fr2.lobe_pa[lobe_cnt] = pa;
                    ++lobe_cnt;
                } else if (structure == 3) {
                    // hotspot
                    assert(hotspot_cnt < 2);
                    fr2.hotspot_flux[hotspot_cnt] = flux;
                    fr2.hotspot_ra[hotspot_cnt] = ra - ra_min;
                    fr2.hotspot_dec[hotspot_cnt] = dec - dec_min;
                    ++hotspot_cnt;

                    if (hotspot_cnt == 2 && simulate_switches[FRII]) {
                        /*
                        * all structures are read in
                        */
                        /* core */
                        int x = (int)round(fr2.core_ra/pix_deg);
                        int y = (int)round(fr2.core_dec/pix_deg);
                        /* hotspot 1 */
                        int x_h1 = (int)round(fr2.hotspot_ra[0]/pix_deg);
                        int y_h1 = (int)round(fr2.hotspot_dec[0]/pix_deg);
                        /* hotspot 2 */
                        int x_h2 = (int)round(fr2.hotspot_ra[1]/pix_deg);
                        int y_h2 = (int)round(fr2.hotspot_dec[1]/pix_deg);
                        /*
                        * check whether fit in the simulation area
                        */
                        if (x < 0 || x >= img_size || y < 0 || y >= img_size)
                            continue;
                        if (x_h1 < 0 || x_h1 >= img_size || y_h1 < 0 || y_h1 >= img_size)
                            continue;
                        if (x_h2 < 0 || x_h2 >= img_size || y_h2 < 0 || y_h2 >= img_size)
                            continue;

                        /*
                        * simulate image
                        */
                        for (vector<double>::size_type k = 0;
                             k < freqs.size();
                             ++k) {
                            /* core */
                            images[k](x, y) += calc_Tb(
                                fr2_core_spec(fr2.core_flux, freqs[k]),
                                pix_area, freqs[k]);
                            /* hotspot 1 */
                            images[k](x_h1, y_h1) += calc_Tb(
                                fr2_hotspot_spec(fr2.hotspot_flux[0], freqs[k]),
                                pix_area, freqs[k]);
                            /* hotspot 2 */
                            images[k](x_h2, y_h2) += calc_Tb(
                                fr2_hotspot_spec(fr2.hotspot_flux[1], freqs[k]),
                                pix_area, freqs[k]);
                        }

                        /* lobe 1 */
                        int xmin, xmax, ymin, ymax;
                        double x_l1 = fr2.lobe_ra[0] / pix_deg;
                        double y_l1 = fr2.lobe_dec[0] / pix_deg;

                        double a1 = 0.5 * fr2.lobe_major_axis[0] / pix_size;
                        double b1 = 0.5 * fr2.lobe_minor_axis[0] / pix_size;
                        double c1 = sqrt(a1*a1 - b1*b1);
                        double s_lobe1 = a1*b1 * M_PI*pix_area;
                        if (s_lobe1 < pix_area)
                            s_lobe1 = pix_area;

                        double f11x = x_l1 + c1*sin(fr2.lobe_pa[0]);
                        double f11y = y_l1 + c1*cos(fr2.lobe_pa[0]);
                        double f21x = x_l1 - c1*sin(fr2.lobe_pa[0]);
                        double f21y = y_l1 - c1*cos(fr2.lobe_pa[0]);

                        xmin = (int) round(x_l1 - a1);
                        xmax = (int) round(x_l1 + a1);
                        ymin = (int) round(y_l1 - a1);
                        ymax = (int) round(y_l1 + a1);
                        for (int i = xmin; i <= xmax; ++i) {
                            for (int j = ymin; j <= ymax; ++j) {
                                if (i < 0 || i >= img_size || j < 0 || j >= img_size)
                                    continue;

                                if (in_ellipse(i, j, f11x, f11y, f21x, f21y, a1)) {
                                    for (vector<double>::size_type k = 0;
                                         k < freqs.size();
                                         ++k) {
                                        images[k](i, j) += calc_Tb(
                                            fr2_lobe_spec(fr2.lobe_flux[0], freqs[k]),
                                            s_lobe1 ,freqs[k]);
                                    }
                                }
                            }
                        }

                        /* lobe 2 */
                        double x_l2 = fr2.lobe_ra[1] / pix_deg;
                        double y_l2 = fr2.lobe_dec[1] / pix_deg;

                        double a2 = 0.5 * fr2.lobe_major_axis[1] / pix_size;
                        double b2 = 0.5 * fr2.lobe_minor_axis[1] / pix_size;
                        double c2 = sqrt(a2*a2 - b2*b2);
                        double s_lobe2 = a2*b2 * M_PI*pix_area;
                        if (s_lobe2 < pix_area)
                            s_lobe2 = pix_area;

                        double f12x = x_l2 + c2*sin(fr2.lobe_pa[1]);
                        double f12y = y_l2 + c2*cos(fr2.lobe_pa[1]);
                        double f22x = x_l2 - c2*sin(fr2.lobe_pa[1]);
                        double f22y = y_l2 - c2*cos(fr2.lobe_pa[1]);

                        xmin = (int) round(x_l2 - a2);
                        xmax = (int) round(x_l2 + a2);
                        ymin = (int) round(y_l2 - a2);
                        ymax = (int) round(y_l2 + a2);
                        for (int i = xmin; i <= xmax; ++i) {
                            for (int j = ymin; j <= ymax; ++j) {
                                if (i < 0 || i >= img_size || j < 0 || j >= img_size)
                                    continue;

                                if (in_ellipse(i, j, f12x, f12y, f22x, f22y, a2)) {
                                    for (vector<double>::size_type k = 0;
                                         k < freqs.size();
                                         ++k) {
                                        images[k](i, j) += calc_Tb(
                                            fr2_lobe_spec(fr2.lobe_flux[1], freqs[k]),
                                            s_lobe2 ,freqs[k]);
                                    }
                                }
                            }
                        }

                        ptr_cnt++;
                        csvout << source_type << "," << redshift << ","
                               << x << "," << y << ","
                               << fr2.core_flux << "," << "," << ","
                               << x_l1 << "," << y_l1 << ","
                               << fr2.lobe_flux[0] << ","
                               << a1*2 << "," << b1*2 << ","
                               << rad2deg(fr2.lobe_pa[0]) << ","  /* lobe 1 */
                               << x_l2 << "," << y_l2 << ","
                               << fr2.lobe_flux[1] << ","
                               << a2*2 << "," << b2*2 << ","
                               << rad2deg(fr2.lobe_pa[1]) << ","  /* lobe 2 */
                               << x_h1 << "," << y_h1 << ","
                               << fr2.hotspot_flux[0] << ","
                               << x_h2 << "," << y_h2 << ","
                               << fr2.hotspot_flux[1]
                               << endl;
                        /* END: simulate image */
                    }
                }
            }

            /*
            * FR I AGN
            */
            if (source_type == FRI) {
                assert(structure==1 || structure==2);
                if (structure == 1) {
                    fr1.reset();
                    fr1.core_ra = ra - ra_min;
                    fr1.core_dec = dec - dec_min;
                    fr1.core_flux = flux;
                } else if (structure == 2) {
                    assert(lobe_cnt < 2);
                    fr1.lobe_ra[lobe_cnt] = ra - ra_min;
                    fr1.lobe_dec[lobe_cnt] = dec - dec_min;
                    fr1.lobe_flux[lobe_cnt] = flux;
                    fr1.lobe_pa[lobe_cnt] = pa;
                    fr1.lobe_major_axis[lobe_cnt] = major_axis;
                    fr1.lobe_minor_axis[lobe_cnt] = minor_axis;
                    ++lobe_cnt;

                    if (lobe_cnt == 2 && simulate_switches[FRI]) {
                        /*
                        * all structures are read in
                        */
                        /* core */
                        int x = (int) round(fr1.core_ra / pix_deg);
                        int y = (int) round(fr1.core_dec / pix_deg);
                        /* lobe 1 */
                        double x_l1 = fr1.lobe_ra[0] / pix_deg;
                        double y_l1 = fr1.lobe_dec[0] / pix_deg;
                        /* lobe 2 */
                        double x_l2 = fr1.lobe_ra[1] / pix_deg;
                        double y_l2 = fr1.lobe_dec[1] / pix_deg;
                        /*
                        * check whether fit in the simulation area
                        */
                        if (x < 0 || x >= img_size || y < 0 || y >= img_size)
                            continue;
                        if (x_l1 < 0 || x_l1 >= img_size || y_l1 < 0 || y_l1 >= img_size)
                            continue;
                        if (x_l2 < 0 || x_l2 >= img_size || y_l2 < 0 || y_l2 >= img_size)
                            continue;

                        /*
                        * simulate image
                        */
                        /* core */
                        for (vector<double>::size_type k = 0;
                             k < freqs.size();
                             ++k) {
                            images[k](x, y) += calc_Tb(
                                fr1_core_spec(fr1.core_flux, freqs[k]),
                                pix_area, freqs[k]);
                        }

                        /* lobe 1 */
                        int xmin, xmax, ymin, ymax;
                        double a1 = 0.5 * fr1.lobe_major_axis[0] / pix_size;
                        double b1 = 0.5 * fr1.lobe_minor_axis[0] / pix_size;
                        double c1 = sqrt(a1*a1 - b1*b1);
                        double s_lobe1 = a1*b1 * M_PI*pix_area;
                        if (s_lobe1 < pix_area)
                            s_lobe1 = pix_area;

                        double f11x = x_l1 + c1*sin(fr1.lobe_pa[0]);
                        double f11y = y_l1 + c1*cos(fr1.lobe_pa[0]);
                        double f21x = x_l1 - c1*sin(fr1.lobe_pa[0]);
                        double f21y = y_l1 - c1*cos(fr1.lobe_pa[0]);

                        xmin = (int) round(x_l1 - a1);
                        xmax = (int) round(x_l1 + a1);
                        ymin = (int) round(y_l1 - a1);
                        ymax = (int) round(y_l1 + a1);
                        for (int i = xmin; i <= xmax; ++i) {
                            for (int j = ymin; j <= ymax; ++j) {
                                if (i < 0 || i >= img_size || j < 0 || j >= img_size)
                                    continue;

                                if (in_ellipse(i, j, f11x, f11y, f21x, f21y, a1)) {
                                    for (vector<double>::size_type k = 0;
                                         k < freqs.size();
                                         ++k) {
                                        images[k](i, j) += calc_Tb(
                                            fr1_lobe_spec(fr1.lobe_flux[0], freqs[k]),
                                            s_lobe1, freqs[k]);
                                    }
                                }
                            }
                        }

                        /* lobe 2 */
                        double a2 = 0.5 * fr1.lobe_major_axis[1] / pix_size;
                        double b2 = 0.5 * fr1.lobe_minor_axis[1] / pix_size;
                        double c2 = sqrt(a2*a2 - b2*b2);
                        double s_lobe2 = a2*b2 * M_PI*pix_area;
                        if (s_lobe2 < pix_area)
                            s_lobe2 = pix_area;

                        double f12x = x_l2 + c2*sin(fr1.lobe_pa[1]);
                        double f12y = y_l2 + c2*cos(fr1.lobe_pa[1]);
                        double f22x = x_l2 - c2*sin(fr1.lobe_pa[1]);
                        double f22y = y_l2 - c2*cos(fr1.lobe_pa[1]);

                        xmin = (int) round(x_l2 - a2);
                        xmax = (int) round(x_l2 + a2);
                        ymin = (int) round(y_l2 - a2);
                        ymax = (int) round(y_l2 + a2);
                        for (int i = xmin; i <= xmax; ++i) {
                            for (int j = ymin; j <= ymax; ++j) {
                                if (i < 0 || i >= img_size || j < 0 || j >= img_size)
                                    continue;

                                if (in_ellipse(i, j, f12x, f12y, f22x, f22y, a2)) {
                                    for (vector<double>::size_type k = 0;
                                         k < freqs.size();
                                         ++k) {
                                        images[k](i, j) += calc_Tb(
                                            fr1_lobe_spec(fr1.lobe_flux[1], freqs[k]),
                                            s_lobe2, freqs[k]);
                                    }
                                }
                            }
                        }

                        ptr_cnt++;
                        csvout << source_type << "," << redshift << ","
                               << x << "," << y << ","
                               << fr1.core_flux << "," << "," << ","
                               << x_l1 << "," << y_l1 << ","
                               << fr1.lobe_flux[0] << ","
                               << a1*2 << "," << b1*2 << ","
                               << rad2deg(fr1.lobe_pa[0]) << ","  /* lobe 1 */
                               << x_l2 << "," << y_l2 << ","
                               << fr1.lobe_flux[1] << ","
                               << a2*2 << "," << b2*2 << ","
                               << rad2deg(fr1.lobe_pa[1])  /* lobe 2 */
                               << endl;
                        /* END: simulate image */
                    }
                }
            }

            old_galaxy = galaxy;
            /* END: one data line */
        }

        ptr_cnt_sum += ptr_cnt;
        /* END: one database file */
    }

    csvout.close();
    cout << "Total simulated point sources: " << ptr_cnt_sum << endl;

    /*
     * Write simulated images to FITS files
     */
    long naxis = 2;
    long naxes[] = {img_size, img_size};
    long nelements = img_size * img_size;
    std::unique_ptr<CCfits::FITS> pfits;

    for (vector<double>::size_type k = 0; k < freqs.size(); ++k) {
        std::ostringstream ss;
        ss << output_prefix
           << std::fixed << std::setprecision(2) << freqs[k]
           << ".fits";
        const std::string fname = ss.str();
        try {
            // "!" to overwrite existing files
            pfits.reset(new CCfits::FITS("!"+fname, DOUBLE_IMG, naxis, naxes));
        }
        catch (CCfits::FITS::CantCreate) {
            cerr << "ERROR: cannot create FITS file: " << fname << endl;
            exit(EXIT_FAILURE);
        }
        long fpixel = 1;
        pfits->pHDU().write(fpixel, nelements, images[k].array());
        pfits->flush();
        cout << "INFO: saved image: " << fname << endl;;
    }
}
