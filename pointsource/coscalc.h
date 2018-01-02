#ifndef CALC_DISTANCE
#define CALC_DISTANCE

extern const double cm;
extern const double s;
extern const double km;
extern const double Mpc;
extern const double kpc;
extern const double yr;
extern const double Gyr;
extern double H;
extern const double c;
//const double c=3e8*100*cm;
extern double omega_m;
extern double omega_l;
extern const double arcsec2arc_ratio;

extern double calc_comoving_distance(double z);
extern double calc_angular_distance(double z);
extern double E(double z);
extern void update_const(double H0,double Om,double Ol);
extern double calc_dvc_dz(double z,double dO);
#endif
//EOF
