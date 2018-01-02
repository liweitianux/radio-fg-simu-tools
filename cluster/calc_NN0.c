/**
 * This tool calculate the synchrotron emission coefficient 'j_nu(t)'
 * at given frequency 'nu' (default: 1.4 GHz) and time 't' (default: 0.0)
 * by 'N0' (default: 1.0) relativistic electron(s) with initial spectral
 * index of Gamma (default: 3.0).
 *
 * The final output column 'NN0':
 *     NN0 = j_{nu=1.4GHz}(t=0) / F_{1.4GHz}(t=0)
 * therefore the flux 'F_{nu}(t)' at different frequency and time can be
 * calculated with:
 *     F_{nu}(t) = j_{nu}(t) / NN0
 *
 * Aaron LI
 * 2015/04/13
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>


/**
 * Custom struct to store the parameters of integrand function
 */
typedef struct {
    double N0;
    double Gamma;
    double B_local;
    double b;
    double age;
    double nu;
} emission_coef_params_t;


/* function prototypes */
double kv( const double x, void *params );
double kv_integral( const double x, const double upper, const double order );
double nu_larmor( const double B_local );
double nu_critical( const double B_local, const double gamma );
double decay_rate( const double B_local, const double z );
double single_electron_spec( const double gamma, const double B_local,
                             const double nu );
double electron_distribution( const double gamma, const double age,
                              const double N0, const double Gamma,
                              const double b );
double emission_coef_integrand( const double gamma, void *params );
double emission_coef( const double gamma_lower, const double gamma_upper,
                      const double N0, const double Gamma,
                      const double B_local, const double b,
                      const double age, const double nu );


/* constants */
const double pi  = 3.14159265;
/* CGS unit system */
const double e   = 4.80320427e-10; /* electric charge (Fr) (ESU, Gaussian) */
const double me  = 9.10938215e-28; /* mass of electron (g) */
const double c   = 2.99792458e10; /* speed of light (cm/s) */
const double yr  = 365 * 24 * 60 * 60; /* 1 yr = ? s */
const double Gyr = 1e9 * 365*24*60*60; /* 1 Gyr = ? s */
const double Gamma = 3.0; /* initial index of the electron energy spectrum */

/* Order of modified Bessel function used here */
const double order = 5.0/3.0;
/* The upper limit of K_{nu} integration. (>700 may cause underflow) */
const double upper = 300.0;

/* Required relative error */
const double epsrel = 1.0e-12;
/* Gauss-Kronrod rule of gsl_integration_qag */
const int gsl_qag_key = 2;

/* The lower & upper limit of emission coefficient integration over gamma */
const double gamma_lower = 1.0;
const double gamma_upper = 1.0e6;


int
main( int argc, char **argv )
{
    if ( argc != 4 ) {
        fprintf( stderr, "Usage:\n" );
        fprintf( stderr, "    %s <B_local_uG> "
              "<in: ra_dec_z_m_F1400_rcDA_T_alpha.dat> "
              "<out: ra_dec_z_m_F1400_rcDA_T_alpha_NN0.dat>\n", argv[0] );
        exit( EXIT_FAILURE );
    }

    double B_local = atof( argv[1] ) * 1e-6; /* uG -> Gauss */

    FILE *fdata = fopen( argv[2], "r" );
    if ( fdata == NULL ) {
        fprintf( stderr, "ERROR: cannot open file %s\n", argv[2] );
        exit( EXIT_FAILURE );
    }
    FILE *fout = fopen( argv[3], "w" );
    if ( fout == NULL ) {
        fprintf( stderr, "ERROR: cannot open file %s\n", argv[3] );
        exit( EXIT_FAILURE );
    }

    double N0 = 1.0; /* normalization coefficient */
    double age = 0.0; /* i.e., to calculate the current state */
    double ra, dec, z, m, F1400, rcDA, T, alpha;
    double b, nu_z, j_nu;

    int i = 0;
    while ( fscanf(fdata, "%lf %lf %lf %lf %lf %lf %lf %lf",
                &ra, &dec, &z, &m, &F1400, &rcDA, &T, &alpha) == 8 ) {
        i++;
        printf("# %d #\n", i);
        b = decay_rate( B_local, z );
        nu_z = 1.4e9 * (1.0 + z); /* 1.4 GHz */
        j_nu = emission_coef( gamma_lower, gamma_upper, N0, Gamma,
                B_local, b, age, nu_z ) / F1400;
        fprintf( fout, "%.10lg %.10lg %.10lg %.10lg %.10lg %.10lg %.10lg %.10lg %.10lg\n",
                ra, dec, z, m, F1400, rcDA, T, alpha, j_nu );
    }

    fclose( fdata );
    fclose( fout );
    return 0;
}


/**
 * Modified Bessel function of the second kind .
 * Irregular modified Bessel function - fractional order (gsl)
 *
 * Parameters:
 *   params: (void *) pointer, pass 'order' through this.
 *           gsl_function require this prototype.
 */
double
kv( const double x, void *params )
{
    double order = *(double *) params;
    return gsl_sf_bessel_Knu( order, x );
}


/**
 * The integration result of K_{nu}, which will be used to
 * calculate the spectrum of one single electron.
 *
 * Parameter:
 *   x: nu / nu_critical
 *   upper: the upper limit of integration.
 *          NOTE: If set a large upper limit, e.g., >700,
 *                exp() underflow may happen.
 *   order: the order of 'K_{nu}' function
 */
double
kv_integral( const double x, const double upper, const double order )
{
    double value, error;
    static double tmp_result;
    static int calculated = 0;

    /* GSL QAG adaptive integration */
    size_t w_size = 100;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc( w_size );
    gsl_function F;
    F.function = kv;
    F.params = (void *) &order;

    if ( x >= 10.0 ) {
        value = sqrt(0.5 * pi * x) * exp(-x);
    } else if ( x >= 0.007 ) {
        gsl_integration_qag( &F, x, upper, 0.0, epsrel,
                w_size, gsl_qag_key, w, &value, &error );
        value *= x;
    } else if ( x > 0.0 ) {
        if ( calculated ) {
            /* Just use the previously save result to save time */
            //fprintf(stderr, "INFO: use previously saved data.\n");
            value = tmp_result * 0.007 * pow(x/0.007, 1.0/3.0);
        } else {
            gsl_integration_qag( &F, 0.007, upper, 0.0, epsrel,
                    w_size, gsl_qag_key, w, &tmp_result, &error );
            value = tmp_result * 0.007 * pow(x/0.007, 1.0/3.0);
            /* Set flag */
            //fprintf(stderr, "INFO: first time calculation; set flag.\n");
            calculated = 1;
        }
    } else {
        value = 0.0;
    }

    gsl_integration_workspace_free( w );
    return value;
}


/**
 * Larmor frequency
 * NOTE: CGS units!
 */
double
nu_larmor( const double B_local )
{
    return e*B_local / (2.0*pi*me*c);
}


/**
 * The critical frequency used in the integration.
 *
 * Parameters:
 *   gamma: Lorentz factor of electron
 */
double
nu_critical( const double B_local, const double gamma )
{
    return 1.5 * nu_larmor(B_local) * pow(gamma, 2);
}


/**
 * Decay rate of the electron spectrum.
 */
double
decay_rate( const double B_local, const double z )
{
    /* Equivalent magnetic field of CMB scattering */
    double B_cmb = pow((1.0+z), 2) * 3.2e-6;
    return 3e-8 * (pow(B_local, 2) + pow(B_cmb, 2)) / (8.0*pi);
}


/**
 * Spectrum of a single electron.
 * XXX: reference??
 *
 * Parameters:
 *   gamma: Lorentz factor of relativistic electron
 *   B_local: local magnetic field (induction)
 *   nu: frequency
 */
double
single_electron_spec( const double gamma, const double B_local,
                      const double nu )
{
    double x = nu / nu_critical(B_local, gamma);
    return 2.0*pi * sqrt(3.0) * pow(e, 2) * nu_larmor(B_local)
        * kv_integral(x, upper, order) / c;
}


/**
 * The number distribution of electron at enery 'gamma' (Lorentz factor)
 * and time 'age'.
 * Calculate the spectral evolution of relativistic electron
 * under synchrotron radiation in the Galactic magnetic field.
 *
 *   N(\gamma, t) d\gamma =
 *       (a) N0 \gamma^{-Gamma} (1-b\gamma t)^{\Gamma-2} d\gamma,
 *       (b) 0, (gamma < gamma_lower, or gamma > gamma_upper)
 *
 * Parameters:
 *   N0: Normalization coefficient
 *   Gamma: initial spectral index of electron energy spectrum
 *   b: electron energy decay rate
 *   gamma: Lorentz factor of electron
 *   age: evolution time (unit: s)
 */
double
electron_distribution( const double gamma, const double age, const double N0,
                       const double Gamma, const double b )
{
    if ( b * gamma * age >= 1.0 ) {
        return 0.0;
    } else {
        return N0 * pow(gamma, -Gamma) * pow((1.0 - b*gamma*age), Gamma-2);
    }
}


/**
 * The integrand function of calculation of emission coefficient.
 */
double
emission_coef_integrand( const double gamma, void *params )
{
    emission_coef_params_t *p = (emission_coef_params_t *) params;
    double N0      = p->N0;
    double Gamma   = p->Gamma;
    double B_local = p->B_local;
    double b       = p->b;
    double age     = p->age;
    double nu      = p->nu;
    return electron_distribution( gamma, age, N0, Gamma, b )
            * single_electron_spec( gamma, B_local, nu );
}


/**
 * Calculate the emission coefficient 'j_{nu}(t)' at frequency 'nu'
 * and time 'age'.
 *
 * Parameters:
 *   gamma_lower, gamma_upper: integration limits
 *   N0: normalization coefficient of electron distribution
 *   Gamma: initial spectral index of electron energy spectrum
 *   B_local: local magnetic field
 *   b: electron energy decay rate
 *   age: evolution time (unit: s)
 *   nu: frequency (unit: Hz)
 */
double
emission_coef( const double gamma_lower, const double gamma_upper,
               const double N0, const double Gamma,
               const double B_local, const double b, const double age,
               const double nu )
{
    emission_coef_params_t params;
    params.N0      = N0;
    params.Gamma   = Gamma;
    params.B_local = B_local;
    params.b       = b;
    params.age     = age;
    params.nu      = nu;

    gsl_function F;
    F.function = emission_coef_integrand;
    F.params = (void *) &params;

    /* GSL QAG adaptive integration */
    size_t w_size = 250;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc( w_size );

    /**
     * Provide points for calculation, otherwise the value of
     * 'b*gamma*age' >= 1.0 when 'gamma' >~ 300.
     */
    size_t i, npts = 31;
    double prop_coef = pow( gamma_upper/gamma_lower, 1.0/(npts-1) );
    double *pts = (double *) malloc( npts * sizeof(double) );
    if ( pts == NULL ) {
        fprintf( stderr, "ERROR: cannot allocate memory for pts.\n" );
        exit( EXIT_FAILURE );
    }
    for ( i=0; i<npts; i++ ) {
        pts[i] = gamma_lower * pow(prop_coef, i);
    }

    gsl_error_handler_t *old_handler;
    /* Turns off the error handler, let program continue after any error. */
    old_handler = gsl_set_error_handler_off();

    int run = 1, ret = 0, j = 0;
    double epsabs = 1.0e-30 * N0;
    double tmp_result = 0.0, value = 0.0, error = 0.0;

    do {
        ret = gsl_integration_qagp( &F, pts, npts, epsabs, 0.0,
                w_size, w, &tmp_result, &error );
        if ( ret == 0 ) {
            value = tmp_result;
            run = 1;
        } else {
            fprintf( stderr,
                    "DEBUG: ret: %d; value: %lg, error: %lg; intervals: %zu\n",
                    ret, tmp_result, error, w->size );
            run = 0; /* Error occurred, stop integration */
            if ( j == 0 ) {
                fprintf( stderr, "WARNING: Error occurred at the "
                        "*first* integration iteration.\n" );
                /* Error occurred at the first run. Retry! */
                value = tmp_result; /* Keep the first integration result */
                run = 0;
            }
        }
        j++; /* Record current integration iteration number */
        epsabs *= 1.0e-4; /* Raise the precision requirement */
    } while ( run && error > tmp_result*0.01 && epsabs > 0 );

    /* Restore the original error handler */
    gsl_set_error_handler( old_handler );
    gsl_integration_workspace_free( w );
    free( pts );
    return value;
}

