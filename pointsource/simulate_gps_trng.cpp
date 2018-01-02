/*-
 * Simulate the GPS/CSS point sources data
 *
 * Wang et al. 2010
 */

/*-
 * ChangeLogs:
 * 2015/05/07:
 *   Use 'TRNG' to generate random numbers
 */

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <utility>

#include <trng/yarn3.hpp>
#include <trng/uniform_dist.hpp>
#include <trng/uniform01_dist.hpp>

#include "coscalc.h"

using namespace std;

const double pi=3.141592654;
const double log10L_min=20;
const double log10L_max=30.5;
const double z_min=0.001;
const double z_max=20;
const double param_v=(z_max-z_min)*(log10L_max-log10L_min);

long nsamples=0L;
double sum_values=0;
double P_max=2.8E2;
double total_num=0;

double calc_dL(double z);
double calc_dL1(double z);
double calc_dA1(double z);
double calc_dVc1(double z);
double dN_dD(double D);
double dN_dlnD(double lnD);
double dr1_dr(double z);
double Lc(double z);
double rho(double log10L,double z);
double rho1(double log10L,double z);
double rho_z(double log10L,double z,double dO);
double rho_z_our(double log10L,double z,double dO);

double get_D(void);
pair<double,double> get_flux_z(void);
double calc_nup_GHz(double Dtrue,double z);


// random number engine
trng::yarn3 R;
// random number distribution
trng::uniform01_dist<> u01;


int main(int argc, char **argv)
{
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <out: gps.dat>" << endl;
        exit(EXIT_FAILURE);
    }

    ofstream ofs(argv[1]);

    R.seed( time(NULL) ); // seed random number engine
    trng::uniform_dist<> ud(-5.0, 5.0);

    /*  
    cout<<rho(25,5)<<endl;
    cout<<rho_z(25,5,.03)<<endl;
    cout<<calc_dVc1(5)<<endl;
    return 0;
    */

    total_num=10;
    for (int i=0; i<10; ++i) {
        get_flux_z();
    }

    pair<double,double> p;
    double F_5GHz, z, Dtrue, da_Mpc, nup_GHz;
    double tk, tn, a1, a2, F_nup;
    double x, y;

    for (int i=0;i<total_num;++i) {
        if ( i%1000==0 ) {
            cerr << i << "\t" << total_num << endl;
        }

        p = get_flux_z();
        F_5GHz = p.first;
        if ( F_5GHz > 1e-6 ) {
            z = p.second;
            da_Mpc = calc_angular_distance(z)/Mpc;
            Dtrue = get_D();
            //double size_arcmin=Dtrue/da_kpc/3.1415*180*60;
            nup_GHz = calc_nup_GHz(Dtrue, z);

            // => Wang et al. 2010, ApJ, Eq. 13 & Eq. 14
            tk = 0.03 * (u01(R) * 2 - 1) + 0.51; // 0.51 +/- 0.03
            tn = 0.06 * (u01(R) * 2 - 1) - 0.73; // -0.73 +/- 0.06
            F_nup = F_5GHz * (1-exp(-1.0)) / pow(5./nup_GHz, tk) /
                    (1-exp(-pow(5./nup_GHz,tn-tk)));

            x = ud(R);
            y = ud(R);
            a1 = 0.05 * (u01(R) * 2 - 1) - 0.21; // -0.21 +/- 0.05
            a2 = 0.05 * (u01(R) * 2 - 1) - 0.65; // -0.65 +/- 0.05
            
            ofs << F_5GHz << "\t" << z << "\t" << "\t"
                << x << "\t" << y << "\t" << Dtrue << "\t"
                << da_Mpc << "\t" << F_nup << "\t "<< nup_GHz << "\t"
                << tk << "\t" << tn << "\t" << a1 << "\t" << a2 << endl;
        }
    }

    return 0;
}


double calc_dL(double z)
{
    //update_const(70*km/s/Mpc,.3,.7);
    return calc_angular_distance(z)*(1+z)*(1+z);
}

double calc_dL1(double z)
{
    const double c=2.99792458E10*cm/s;
    const double H1=50*km/s/Mpc;
    return c*(1+z)/H1*2*(1-1/sqrt(1+z));
}

double calc_dA1(double z)
{
    return calc_dL1(z)/(1+z)/(1+z)/Mpc;
}

double calc_dVc1(double z)
{
    const double c=2.99792458E5;
    const double H1=50;
    double E=sqrt((1+z)*(1+z)*(1+z));
    double da1=calc_dA1(z);
    return c/H1*(1+z)*(1+z)*da1*da1/E;
}

double dN_dD(double D)
{
    if (D<=.5&&D>=.01) {
        return 1/D;
    } else if (D>.5&&D<=20) {
        return pow(D,-.6)/pow(.5,-.6)*1/.5;
    }
    return 0;
}

double dN_dlnD(double lnD)
{
    double D=exp(lnD);
    return D*dN_dD(D);
}

double dr1_dr(double z)
{
    const double H0=71;
    const double H1=50;
    const double Om=.27;
    const double Ol=.73;
    return H0*sqrt(Om*(1+z)*(1+z)*(1+z)+Ol) / (H1*pow(1+z,1.5));
}

double Lc(double z)
{
    const double a0=25.99;
    const double a1=1.26;
    const double a2=-.26;
    return pow(10.,a0+a1*z+a2*z*z);
}

double rho(double log10L,double z)
{
    return rho1(log10L,z)*pow(calc_dL1(z)/calc_dL(z),2)*dr1_dr(z);
}

double rho1(double log10L,double z)
{
    const double rho0=pow(10.,-6.91);
    const double a=.69;
    const double b=2.17;
    double L1=pow(10.,log10L);
    return rho0/(pow(4*pi*L1/Lc(z),a)+pow(4*pi*L1/Lc(z),b));
}

double rho_z(double log10L,double z,double dO=1)
{
    //  update_const(50*km/s/Mpc,.0,.0);
    return rho(log10L,z)*calc_dVc1(z)*dO;
}

double rho_z_our(double log10L,double z,double dO=1)
{
    double dL1_dL=calc_dL1(z)/calc_dL(z);
    double log10L1=log10L+2*log10(dL1_dL);
    return rho_z(log10L1,z,dO);
}

double get_D(void)
{
    const double lnD_min=log(.01);
    const double lnD_max=log(20);
    const double P_max=dN_dlnD(lnD_min);
    double lnD;
    do {
        lnD=u01(R)*(lnD_max-lnD_min)+lnD_min;
        //cout<<dN_dlnD(lnD)<<"\t"<<P_max<<endl;
    } while (dN_dlnD(lnD)<P_max*u01(R));
  
    return exp(lnD);
}

pair<double,double> get_flux_z(void)
{
    double z,log10L,P;
    double P1=0;
    bool reject=true;
    do {
        reject=false;
        z=u01(R)*(z_max-z_min)+z_min;
        log10L=u01(R)*(log10L_max-log10L_min)+log10L_min;
        P=rho_z_our(log10L,z);
        P1=P_max*u01(R);
        //cout<<P<<"\t"<<P1<<"\t"<<endl;
        ++nsamples;
        sum_values+=P;
        total_num=sum_values/nsamples*param_v*0.03046;
        if (P>=P_max) {
            reject=true;
            cerr<<"updated:"<<P<<"\t"<<P_max<<endl;
            P_max=P;
        }
    } while (P<P1||reject);
    double L=pow(10.,log10L);
    double dL=calc_dL(z);
    double flux=L/(dL/100)/(dL/100)/1E-26;
    return make_pair(flux,z);
}

double calc_nup_GHz(double Dtrue, double z)
{
    // log( nup (1+z) / GHz) = -0.21 (+/- 0.05) - 0.65 (+/- 0.05) log(D / kpc)
    return pow(10.,-.21-.65*log10(Dtrue))/(1+z);
}

