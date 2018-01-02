#ifndef F_ROSAT_H
#define F_ROSAT_H

extern double f_bol_Tz(double T,double z);
extern double L_rosat_Tz(double T,double z);
extern double f_rosat_Tz(double T,double z);
extern double f_rosat_z(double z);
extern double get_bolo_flux(double z,double T);
extern double get_soft_flux(double z,double T);
extern double get_s2b_rate(double z,double T);
extern double current_T;

#endif
//EOF



