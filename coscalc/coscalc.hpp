#ifndef COSCALC_HPP
#define COSCALC_HPP

static const int CGS=1;
static const int SI=2;

extern double cm;
extern double gram;
extern double sec;
extern double erg;

extern double meter;
extern double kg;
extern double km;


extern double hr;
extern double day;
extern double yr;
extern double kpc;
extern double Mpc;
extern double J;


extern double mp;
extern double me;
extern double eV;
extern double keV;

extern double H0;



extern double E(double z,double omega_m=.27,double omega_l=.73);
extern double calc_dc(double z,double omega_m=.27,double omega_l=.73);
extern double calc_dm(double z,double omega_m=.27,double omega_l=.73);
extern double calc_da(double z,double omega_m=.27,double omega_l=.73);
extern double calc_dl(double z,double omega_m=.27,double omega_l=.73);
#endif
//EOF
