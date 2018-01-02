/*-
 * Simulate the GPS/CSS point sources data
 *
 * Wang et al. 2010
 */

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <fstream>
#include <utility>

#include "coscalc.h"

using namespace std;

const double pi=3.141592653;
const double log10L_min=20;
const double log10L_max=30.5;
const double z_min=0.001;
const double z_max=20;
const double param_v=(z_max-z_min)*(log10L_max-log10L_min);

size_t nsamples=0;
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
double random_norm();
double random_uniform(double a,double b);

double get_D(void);
pair<double,double> get_flux_z(void);
double calc_vp_GHz(double Dtrue,double z);
double calc_gps_spec(double s_5000,double vp_GHz,double nu_MHz);


int main(int argc, char **argv)
{
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <out: gps.dat>" << endl;
        exit(EXIT_FAILURE);
    }

    ofstream ofs(argv[1]);

    srand(time(NULL));
    total_num=10;

    /*  
    cout<<rho(25,5)<<endl;
    cout<<rho_z(25,5,.03)<<endl;
    cout<<calc_dVc1(5)<<endl;
    return 0;
    */

    for (int i=0; i<10; ++i) {
        get_flux_z();
    }

    pair<double,double> p;
    for (int i=0;i<total_num;++i) {
        if (i%1000==0) {
            cerr << i << "\t" << total_num << endl;
        }

        p=get_flux_z();
        if (p.first>1E-6) {
            double flux_5000=p.first;
            double z=p.second;
            double Dtrue=get_D();
            double da_Mpc=calc_angular_distance(z)/Mpc;
            //double size_arcmin=Dtrue/da_kpc/3.1415*180*60;
            double vp_GHz=calc_vp_GHz(Dtrue,z);

            double k=((double)rand())/RAND_MAX*0.03*2-0.03+0.51;
            double l=((double)rand())/RAND_MAX*0.06*2-0.06-0.73;
            double a1=((double)rand())/RAND_MAX*0.05*2-0.05-0.21;
            double a2=((double)rand())/RAND_MAX*0.05*2-0.05-0.65;
            
            double S_vp=flux_5000/(1/(1-exp(-1.))*pow(5./vp_GHz,k)*(1-exp(-pow(5./vp_GHz,l-k))));
            
            double x=random_uniform(-5.,5.);
            double y=random_uniform(-5.,5.);
            
            ofs << flux_5000 << "\t" << z << "\t" << "\t"
                << x << "\t" << y << "\t" << Dtrue << "\t"
                << da_Mpc << "\t" << S_vp << "\t "<< vp_GHz << "\t"
                << k << "\t" << l << "\t" << a1 << "\t" << a2 << endl;
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

double random_norm()
{
    return (double) rand() / RAND_MAX;
}

double random_uniform(double a,double b)
{
    return random_norm()*(b-a)+a;
}

double get_D(void)
{
    const double lnD_min=log(.01);
    const double lnD_max=log(20);
    const double P_max=dN_dlnD(lnD_min);
    double lnD;
    do {
        lnD=random_norm()*(lnD_max-lnD_min)+lnD_min;
        //cout<<dN_dlnD(lnD)<<"\t"<<P_max<<endl;
    } while (dN_dlnD(lnD)<P_max*random_norm());
  
    return exp(lnD);
}

pair<double,double> get_flux_z(void)
{
    double z,log10L,P;
    double P1=0;
    bool reject=true;
    do {
        reject=false;
        z=random_norm()*(z_max-z_min)+z_min;
        log10L=random_norm()*(log10L_max-log10L_min)+log10L_min;
        P=rho_z_our(log10L,z);
        P1=P_max*random_norm();
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

double calc_vp_GHz(double Dtrue,double z)
{
    double vp=pow(10.,-.21-.65*log10(Dtrue))/(1+z);
    return vp;
}

double calc_gps_spec(double s_5000,double vp_GHz,double nu_MHz)
{
    double k=((double)rand())/RAND_MAX*0.03*2-0.03+0.51;
    double l=((double)rand())/RAND_MAX*0.06*2-0.06-0.73;
    //double a1=((double)rand())/RAND_MAX*0.05*2-0.05-0.21;
    //double a2=((double)rand())/RAND_MAX*0.05*2-0.05-0.65;

    double S_vp=s_5000/(1/(1-exp(-1.))*pow(5./vp_GHz,k)*(1-exp(-pow(5./vp_GHz,l-k))));
    double v_GHz=nu_MHz/1000.;
    double Sv=S_vp*(1/(1-exp(-1))*pow(v_GHz/vp_GHz,k)*(1-exp(-pow(v_GHz/vp_GHz,l-k))));
    return Sv;
}

