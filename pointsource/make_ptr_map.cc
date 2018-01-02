/*
 * Simulate maps of point sources based on the database by Wilman et al. 2008
 * http://s-cubed.physics.ox.ac.uk/s3_sex/database_structure
 *
 * Junhua GU, Jingying Wang, 2010
 *
 * Weitian LI, 2017
 */

#include <fstream>
#include <iostream>
#include <cassert>
#include <cmath>
#include <string>
#include <vector>

#include <unistd.h>

#include "fio/include/fio.h"

using namespace std;
using namespace blitz;

/*
 * NOTE: The simulating sky size is fixed to be 10x10 [deg^2], which is
 *       constrained by the obtained database from Wilman et al. 2008.
 *       The obtained database covers R.A. -5:5 [deg] and Dec. -5:5 [deg].
 */
const double sky_size_degree = 10.0;  // Fixed to 10x10 deg^2
const double sky_size_arcsec = sky_size_degree * 3600.0;

const int FRI=1;
const int FRII=2;
const int RQQ=3;
const int SF=4;
const int SB=5;

//int simulate_switches[6]={-1,1,0,0,0,0}; //FRI
//int simulate_switches[6]={-1,0,1,0,0,0}; //FRII
//int simulate_switches[6]={-1,0,0,1,0,0}; //RQQ
//int simulate_switches[6]={-1,0,0,0,1,0}; //SF
//int simulate_switches[6]={-1,0,0,0,0,1}; //SB
int simulate_switches[6]={-1,1,1,1,1,1}; //ALL


inline int mod(int a, int b)
{
  int result;
  result = a % b;
  result = result<0 ? result+b : result;
  return result;
}


struct FRI_info
{
  double core_ra;
  double core_dec;
  double core_flux;
  double lobe_flux[2];
  double lobe_major_axis[2];
  double lobe_minor_axis[2];
  double lobe_ra[2];
  double lobe_dec[2];
  double lobe_position_angle[2];
  
  double redshift;
  double cosv;

  void reset()
  {
    core_ra=-1;
    core_dec=-1;
    core_flux=-1;
    for(int i=0;i<2;++i)
      {
	lobe_flux[i]=-1;
	lobe_major_axis[i]=-1;
	lobe_minor_axis[i]=-1;
	
	lobe_ra[i]=-1;
	lobe_dec[i]=-1;
	lobe_position_angle[i]=-1;
	
      }
    redshift=-1;
    cosv=-1;

  }
};

struct FRII_info
{
  double core_ra;
  double core_dec;
  double core_flux;
  double lobe_flux[2];
  double lobe_major_axis[2];
  double lobe_minor_axis[2];
  double lobe_ra[2];
  double lobe_dec[2];
  double lobe_position_angle[2];
  
  double hotspot_flux[2];
  double hotspot_ra[2];
  double hotspot_dec[2];
  double cosv;
  double redshift;

  void reset()
  {
    core_ra=-1;
    core_dec=-1;
    core_flux=-1;
    for(int i=0;i<2;++i)
      {
	lobe_flux[i]=-1;
	lobe_major_axis[i]=-1;
	lobe_minor_axis[i]=-1;
	lobe_ra[i]=-1;
	lobe_dec[i]=-1;
	lobe_position_angle[i]=-1;
	
	hotspot_flux[i]=-1;
	hotspot_ra[i]=-1;
	hotspot_dec[i]=-1;
	cosv=-1;
	redshift=-1;
      }
 
  }
};

struct RQQ_info
{
  double flux;
  double ra;
  double dec;
  double redshift;
  void reset()
  {
    flux=-1;
    ra=-1;
    dec=-1;
    redshift=-1;
  }
};

struct SF_info
{
  double flux;
  double ra;
  double dec;
  double major_axis;
  double minor_axis;
  double position_angle;
  double m_hi;
  int type;
  double redshift;

  void reset()
  {
    flux=-1;
    ra=-1;
    dec=-1;
    major_axis=-1;
    minor_axis=-1;
    position_angle=-1;
    m_hi=-1;
    redshift=-1;

  }
};


double rqq_spec(double I_151,double freq)
{
  return pow(freq/151E6,-.7)*I_151;
}

double fr1_lobe_spec(double I_151,double freq)
{
  return pow(freq/151E6,-.75)*I_151;
}

double fr1_core_spec(double I_151,double freq)
{
  double a0=log10(I_151)-.07*log10(151E6)+.29*log10(151E6)*log10(151E6);
  double lgS=a0+.07*log10(freq)-.29*log10(freq)*log10(freq);
  return pow(10.,lgS);
}

double fr2_core_spec(double I_151,double freq)
{
  double a0=log10(I_151)-.07*log10(151E6)+.29*log10(151E6)*log10(151E6);
  double lgS=a0+.07*log10(freq)-.29*log10(freq)*log10(freq);
  return pow(10.,lgS);
}

double fr2_lobe_spec(double I_151,double freq)
{
  return pow(freq/151E6,-.75)*I_151;
}

double fr2_hotspot_spec(double I_151,double freq)
{
  return pow(freq/151E6,-.75)*I_151;
}

double sf_spec(double I_151,double freq)
{
  return pow(freq/151E6,-.7)*I_151;
}

double sb_spec(double I_151,double freq)
{
  return pow(freq/151E6,-.7)*I_151;
}

double calc_Tb(double flux_in_Jy,double pixel_area/*in arcsec^2*/,double freq)
{
  const double c=2.99792458E8;
  const double kb=1.38E-23;
  double Omegab=pixel_area/(3600.*180./3.1415)/(3600.*180./3.1415);
  
  double Sb=flux_in_Jy*1E-26/Omegab;
  double result=Sb/2/freq/freq*c*c/kb;
  return result;
}


void usage(const char *name) {
  fprintf(stderr, "Usage: %s [ -o output_prefix ] [ -s image_size ] "
      "-i <dbfiles.list> freqMHz ...\n", name);
  fprintf(stderr, "\n");
  fprintf(stderr, "    -o : output prefix; default: 'ptr_'\n");
  fprintf(stderr, "    -s : image size for the 10x10 deg^2 sky; default: 1800\n");
  fprintf(stderr, "    -i : filename of list of input database files; [required]\n");
  fprintf(stderr, "    freqMHz : list of simulating frequencies [MHz]; [required]\n");
  fprintf(stderr, "\n");
  exit(EXIT_FAILURE);
}


int main(int argc, char *argv[])
{
  const char *progname = argv[0];
  int opt;
  int img_size = 1800;
  const char *output_prefix = "ptr_";
  const char *dbfiles_list = NULL;

  while ((opt = getopt(argc, argv, "i:o:s:")) != -1) {
    switch (opt) {
      case 'i':
        dbfiles_list = optarg;
        break;
      case 'o':
        output_prefix = optarg;
        break;
      case 's':
        img_size = atoi(optarg);
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
  if (dbfiles_list == NULL) {
    usage(progname);
  }

  ifstream input_list(dbfiles_list);
  printf("INFO: database files list: %s\n", dbfiles_list);
  const double pixel_size = sky_size_arcsec / img_size;  // [arcsec]
  const double pixel_area = pixel_size * pixel_size;  // [arcsec^2]
  printf("INFO: image size: %d\n", img_size);
  printf("INFO: pixel size: %.1f [arcsec]\n", pixel_size);
  printf("INFO: sky size: 10x10 [deg^2] (FIXED)\n");

  vector<double> freq_vec;
  vector<blitz::Array<double,2> > img_vec;
  printf("INFO: frequencies [MHz]:");
  for(int i=0; i<argc; ++i)
    {
      freq_vec.push_back(atof(argv[i]));
      img_vec.push_back(blitz::Array<double,2>(img_size,img_size));
      printf(" %.2f", freq_vec.back());
    }
  printf("\n");
  printf("INFO: number of frequencies: %lu\n", freq_vec.size());

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
  double ra,                /* 7. */
         dec,               /* 8. */
         distance,          /* 9. */
         redshift,          /* 10. */
         position_angle,    /* 11. */
         major_axis,        /* 12. */
         minor_axis,        /* 13. */
         i_151,             /* 14. */
         i_610,             /* 15. */
         i_1400,            /* 16. */
         i_4860,            /* 17. */
         i_18000,           /* 18. */
         m_hi,              /* 19. */
         cos_va;            /* 20. */

  for(;;)
    {
      getline(input_list, fname);
      if(!input_list.good())
	{
	  break;
	}
      cout << "Database: " << fname << endl;
      ifstream ifs(fname.c_str());
  
      old_galaxy = -1;
      lobe_cnt = 0;
      hotspot_cnt = 0;

      FRI_info fr1;
      FRII_info fr2;
      SF_info sf;
      RQQ_info rqq;

      for(int n=0;;++n)
	{
	  ifs >> source >> cluster >> galaxy >> sftype >> agntype >> structure;

	  if(!ifs.good())
	    {
	      cout << "Source/component counts: "<< n << endl;
	      break;
	    }

	  ifs >> ra >> dec >> distance >> redshift >> position_angle
	      >> major_axis >> minor_axis >> i_151 >> i_610 >> i_1400
	      >> i_4860 >> i_18000;
 
	  if(sftype)
	    {
	      ifs >> m_hi;
	    }
	  ifs >> cos_va;

	  if(sftype==1)
	    {
	      source_type=SF;
	    }
	  else if(sftype==2)
	    {
	      source_type=SB;
	    }
	  else if(agntype==1)
	    {
	      source_type=RQQ;
	    }
	  else if(agntype==2)
	    {
	      source_type=FRI;
	    }
	  else if(agntype==3)
	    {
	      source_type=FRII;
	    }
	  else
	    {
	      cout<<sftype<<"\t"<<agntype<<endl;
	      assert(0);
	    }

	  if(old_galaxy!=galaxy)  //different component in one galaxy
	    {
	      hotspot_cnt=0;
	      lobe_cnt=0;
	      fr1.reset();
	      fr2.reset();
	      rqq.reset();
	      sf.reset();
	    }

          /*
           * Radio-quiet AGN
           */
	  if(source_type==RQQ)
	    {
	      if(simulate_switches[RQQ])
		{
		  rqq.reset();
		  rqq.ra=ra;
		  rqq.dec=dec;
		  rqq.flux=pow(10.,i_151);
		  rqq.redshift=redshift;
		  int x=(int)(ra/sky_size_degree*img_size);
		  int y=(int)(dec/sky_size_degree*img_size);
		  for(vector<double>::size_type k=0;k<freq_vec.size();++k)
		    {
		      double freq=freq_vec[k]*1E6;
		      img_vec[k](mod(x,img_size),mod(y,img_size))+=calc_Tb(rqq_spec(rqq.flux,freq),pixel_area,freq);
		    }
		}
	    }

          /*
           * Star-formation / starburst galaxies
           */
	  if(source_type==SF||source_type==SB)
	    {
	      sf.reset();
	      assert(structure==4);
	      sf.m_hi=m_hi;
	      sf.ra=ra;
	      sf.dec=dec;
	      sf.flux=pow(10.,i_151);
	      sf.major_axis=major_axis;
	      sf.minor_axis=minor_axis;
	      sf.position_angle=position_angle;
	      sf.redshift=redshift;
	      sf.type=sftype;

	      if((simulate_switches[SF]&&source_type==SF)||
		 (simulate_switches[SB]&&source_type==SB))
		{
		  double a=sf.major_axis/2/sky_size_arcsec*img_size;
		  double b=sf.minor_axis/2/sky_size_arcsec*img_size;
		  
		  double c=sqrt(a*a-b*b);
		  
		  double center_x=(ra/sky_size_degree*img_size);
		  double center_y=(dec/sky_size_degree*img_size);
		  
		  double f1x=center_x+c*sin(position_angle);
		  double f1y=center_y+c*cos(position_angle);
		  
		  double f2x=center_x-c*sin(position_angle);
		  double f2y=center_y-c*cos(position_angle);
		  
		  double s=3.14*major_axis/2.*minor_axis/2.;  // [arcsec^2]
		  if (s < pixel_area)
		    s = pixel_area;
		  
		  //cout<<f1x<<"\t"<<f1y<<"\t"<<center_x<<"\t"<<center_y<<"\t"<<f2x<<"\t"<<f2y<<endl;
		  
		  int xmin = (int)(center_x - a - 1);
		  int xmax = (int)(center_x + a + 1);
		  int ymin = (int)(center_y - a - 1);
		  int ymax = (int)(center_y + a + 1);
		  for(int i = xmin; i <= xmax; ++i)
		    {
		      for(int j = ymin; j <= ymax; ++j)
			{
			  if(sqrt((i-f1x)*(i-f1x)+(j-f1y)*(j-f1y))
			     +sqrt((i-f2x)*(i-f2x)+(j-f2y)*(j-f2y))<2*a)
			    {
			      assert(sf.type);
			      if(sf.type==1)
				{
				  for(vector<double>::size_type k=0;k<freq_vec.size();++k)
				    {
				      double freq=freq_vec[k]*1E6;
				      img_vec[k](mod(i,img_size),mod(j,img_size))+=calc_Tb(sf_spec(sf.flux,freq),s,freq);
				    }
				}
			      else if(sf.type==2)
				{
				   for(vector<double>::size_type k=0;k<freq_vec.size();++k)
				    {
				      double freq=freq_vec[k]*1E6;
				      img_vec[k](mod(i,img_size),mod(j,img_size))+=calc_Tb(sb_spec(sf.flux,freq),s,freq);
				    }
				}
			      else
				{
				  assert(0);
				}
			    }
			}
		    }
		}
	      
	    }

          /*
           * FR II AGN
           */
	  if (source_type==FRII)
	    {
	      assert(structure>=1&&structure<=3);
	      //cout<<structure<<endl;
	      if(structure==1)//core
		{
		  fr1.reset();
		  fr2.core_ra=ra;
		  fr2.core_dec=dec;
		  fr2.core_flux=pow(10.,i_151);
		  fr2.redshift=redshift;
		  fr2.cosv=cos_va;
		}
	      else if(structure==2)//lobe
		{
		  assert(lobe_cnt<2);
		  fr2.lobe_flux[lobe_cnt]=pow(10.,i_151);
		  fr2.lobe_ra[lobe_cnt]=ra;
		  fr2.lobe_dec[lobe_cnt]=dec;
		  fr2.lobe_major_axis[lobe_cnt]=major_axis;
		  fr2.lobe_minor_axis[lobe_cnt]=minor_axis;
		  fr2.lobe_position_angle[lobe_cnt]=position_angle;
		  ++lobe_cnt;
		}
	      else if(structure==3)//hotspot
		{
		  assert(hotspot_cnt<2);
		  fr2.hotspot_flux[hotspot_cnt]=pow(10.,i_151);
		  fr2.hotspot_ra[hotspot_cnt]=ra;
		  fr2.hotspot_dec[hotspot_cnt]=dec;
		  ++hotspot_cnt;
		  if(hotspot_cnt==2&&simulate_switches[FRII])
		    {
		      //output
		      int x_core=(int)(fr2.core_ra/sky_size_degree*img_size);
		      int y_core=(int)(fr2.core_dec/sky_size_degree*img_size);
		      for(vector<double>::size_type k=0;k<freq_vec.size();++k)
			{
			  double freq=freq_vec[k]*1E6;
			  img_vec[k](mod(x_core,img_size),mod(y_core,img_size))+=calc_Tb(fr2_core_spec(fr2.core_flux,freq),pixel_area,freq);
			}
		  
		      int x_hotspot_1=(int)(fr2.hotspot_ra[0]/sky_size_degree*img_size);
		      int y_hotspot_1=int(fr2.hotspot_dec[0]/sky_size_degree*img_size);
		      for(vector<double>::size_type k=0;k<freq_vec.size();++k)
			{
			  double freq=freq_vec[k]*1E6;
			  img_vec[k](mod(x_hotspot_1,img_size),mod(y_hotspot_1,img_size))+=calc_Tb(fr2_hotspot_spec(fr2.hotspot_flux[0],freq),pixel_area,freq);
			}

		      int x_hotspot_2=(int)(fr2.hotspot_ra[1]/sky_size_degree*img_size);
		      int y_hotspot_2=(int)(fr2.hotspot_dec[1]/sky_size_degree*img_size);
		      for(vector<double>::size_type k=0;k<freq_vec.size();++k)
			{
			  double freq=freq_vec[k]*1E6;
			  img_vec[k](mod(x_hotspot_2,img_size),mod(y_hotspot_2,img_size))+=calc_Tb(fr2_hotspot_spec(fr2.hotspot_flux[1],freq),pixel_area,freq);
			}

		      double x_lobe_1=(fr2.lobe_ra[0]/sky_size_degree*img_size);
		      double y_lobe_1=(fr2.lobe_dec[0]/sky_size_degree*img_size);

		      double x_lobe_2=(fr2.lobe_ra[1]/sky_size_degree*img_size);
		      double y_lobe_2=(fr2.lobe_dec[1]/sky_size_degree*img_size);
		  
		      double a1=fr2.lobe_major_axis[0]/2;
		      double a2=fr2.lobe_major_axis[1]/2;
		      double b1=fr2.lobe_minor_axis[0]/2;
		      double b2=fr2.lobe_minor_axis[1]/2;
		      double c1=sqrt(a1*a1-b1*b1);
		      double c2=sqrt(a2*a2-b2*b2);
		      double s_lobe1=3.14*a1*b1;
		      double s_lobe2=3.14*a2*b2;
		      if (s_lobe1 < pixel_area)
			s_lobe1 = pixel_area;
		      if (s_lobe2 < pixel_area)
			s_lobe2 = pixel_area;

		      double f11_x=x_lobe_1+c1*sin(fr2.lobe_position_angle[0])/sky_size_arcsec*img_size;
		      double f11_y=y_lobe_1+c1*cos(fr2.lobe_position_angle[0])/sky_size_arcsec*img_size;

		      double f21_x=x_lobe_1-c1*sin(fr2.lobe_position_angle[0])/sky_size_arcsec*img_size;
		      double f21_y=y_lobe_1-c1*cos(fr2.lobe_position_angle[0])/sky_size_arcsec*img_size;

		      double f12_x=x_lobe_2+c2*sin(fr2.lobe_position_angle[1])/sky_size_arcsec*img_size;
		      double f12_y=y_lobe_2+c2*cos(fr2.lobe_position_angle[1])/sky_size_arcsec*img_size;

		      double f22_x=x_lobe_2-c2*sin(fr2.lobe_position_angle[1])/sky_size_arcsec*img_size;
		      double f22_y=y_lobe_2-c2*cos(fr2.lobe_position_angle[1])/sky_size_arcsec*img_size;

		      //		  cout<<f11_x<<"\t"<<f11_y<<"\t"<<f21_x<<"\t"<<f21_y<<"\t"<<x_lobe_1<<"\t"<<y_lobe_1<<endl;

		      int xmin = (int)(x_lobe_1 - a1 / sky_size_arcsec * img_size);
		      int xmax = (int)(x_lobe_1 + a1 / sky_size_arcsec * img_size);
		      int ymin = (int)(y_lobe_1 - a1 / sky_size_arcsec * img_size);
		      int ymax = (int)(y_lobe_1 + a1 / sky_size_arcsec * img_size);
		      for(int i = xmin; i <= xmax; ++i)
			{
			  for(int j = ymin; j <= ymax; ++j)
			    {
			      if(sqrt((i-f11_x)*(i-f11_x)+(j-f11_y)*(j-f11_y))+
				 sqrt((i-f21_x)*(i-f21_x)+(j-f21_y)*(j-f21_y))<2*a1/sky_size_arcsec*img_size)
				{
				  for(vector<double>::size_type k=0;k<freq_vec.size();++k)
				    {
				      double freq=freq_vec[k]*1E6;
				      img_vec[k](mod(i,img_size),mod(j,img_size))+=calc_Tb(fr2_lobe_spec(fr2.lobe_flux[0],freq),s_lobe1,freq);
				    }
				}
			    }
			}

		      xmin = (int)(x_lobe_2 - a2 / sky_size_arcsec * img_size);
		      xmax = (int)(x_lobe_2 + a2 / sky_size_arcsec * img_size);
		      ymin = (int)(y_lobe_2 - a2 / sky_size_arcsec * img_size);
		      ymax = (int)(y_lobe_2 + a2 / sky_size_arcsec * img_size);
		      for(int i = xmin; i <= xmax; ++i)
			{
			  for(int j = ymin; j <= ymax; ++j)
			    {
			      if(sqrt((i-f12_x)*(i-f12_x)+(j-f12_y)*(j-f12_y))+
				 sqrt((i-f22_x)*(i-f22_x)+(j-f22_y)*(j-f22_y))<2*a2/sky_size_arcsec*img_size)
				{
				  for(vector<double>::size_type k=0;k<freq_vec.size();++k)
				    {
				      double freq=freq_vec[k]*1E6;
				      img_vec[k](mod(i,img_size),mod(j,img_size))+=calc_Tb(fr2_lobe_spec(fr2.lobe_flux[1],freq),s_lobe2,freq);
				    }
				}
			    }
			}

		    }
		}
	    }

          /*
           * FR I AGN
           */
	  if(source_type==FRI)
	    {
	      assert(structure==1||structure==2);
	      if(structure==1)
		{
		  fr1.reset();
		  fr1.core_ra=ra;
		  fr1.core_dec=dec;
		  fr1.core_flux=pow(10.,i_151);
		}
	      else if(structure==2)
		{
		  assert(lobe_cnt<2);
		  fr1.lobe_ra[lobe_cnt]=ra;
		  fr1.lobe_dec[lobe_cnt]=dec;
		  fr1.lobe_flux[lobe_cnt]=pow(10.,i_151);
		  fr1.lobe_position_angle[lobe_cnt]=position_angle;
		  fr1.lobe_major_axis[lobe_cnt]=major_axis;
		  fr1.lobe_minor_axis[lobe_cnt]=minor_axis;
		  ++lobe_cnt;
		  if(lobe_cnt==2&&simulate_switches[FRI])
		    {
		      int x_core=(int)(fr1.core_ra/sky_size_degree*img_size);
		      int y_core=(int)(fr1.core_dec/sky_size_degree*img_size);
		      for(vector<double>::size_type k=0;k<freq_vec.size();++k)
			{
			  double freq=freq_vec[k]*1E6;
			  img_vec[k](mod(x_core,img_size),mod(y_core,img_size))+=calc_Tb(fr1_core_spec(fr1.core_flux,freq),pixel_area,freq);
			}
		  
		      double x_lobe_1=(fr1.lobe_ra[0]/sky_size_degree*img_size);
		      double y_lobe_1=(fr1.lobe_dec[0]/sky_size_degree*img_size);

		      double x_lobe_2=(fr1.lobe_ra[1]/sky_size_degree*img_size);
		      double y_lobe_2=(fr1.lobe_dec[1]/sky_size_degree*img_size);
		  
		      double a1=fr1.lobe_major_axis[0]/2;
		      double a2=fr1.lobe_major_axis[1]/2;
		      double b1=fr1.lobe_minor_axis[0]/2;
		      double b2=fr1.lobe_minor_axis[1]/2;
		      double c1=sqrt(a1*a1-b1*b1);
		      double c2=sqrt(a2*a2-b2*b2);
		      double s_lobe1=3.14*a1*b1;
		      double s_lobe2=3.14*a2*b2;
		      if (s_lobe1 < pixel_area)
			s_lobe1 = pixel_area;
		      if (s_lobe2 < pixel_area)
			s_lobe2 = pixel_area;

		      double f11_x=x_lobe_1+c1*sin(fr1.lobe_position_angle[0])/sky_size_arcsec*img_size;
		      double f11_y=y_lobe_1+c1*cos(fr1.lobe_position_angle[0])/sky_size_arcsec*img_size;

		      double f21_x=x_lobe_1-c1*sin(fr1.lobe_position_angle[0])/sky_size_arcsec*img_size;
		      double f21_y=y_lobe_1-c1*cos(fr1.lobe_position_angle[0])/sky_size_arcsec*img_size;

		      double f12_x=x_lobe_2+c2*sin(fr1.lobe_position_angle[1])/sky_size_arcsec*img_size;
		      double f12_y=y_lobe_2+c2*cos(fr1.lobe_position_angle[1])/sky_size_arcsec*img_size;

		      double f22_x=x_lobe_2-c2*sin(fr1.lobe_position_angle[1])/sky_size_arcsec*img_size;
		      double f22_y=y_lobe_2-c2*cos(fr1.lobe_position_angle[1])/sky_size_arcsec*img_size;

		      //		  cout<<f11_x<<"\t"<<f11_y<<"\t"<<f21_x<<"\t"<<f21_y<<"\t"<<x_lobe_1<<"\t"<<y_lobe_1<<endl;

		      int xmin = (int)(x_lobe_1 - a1 / sky_size_arcsec * img_size);
		      int xmax = (int)(x_lobe_1 + a1 / sky_size_arcsec * img_size);
		      int ymin = (int)(y_lobe_1 - a1 / sky_size_arcsec * img_size);
		      int ymax = (int)(y_lobe_1 + a1 / sky_size_arcsec * img_size);
		      for(int i = xmin; i <= xmax; ++i)
			{
			  for(int j = ymin; j <= ymax; ++j)
			    {
			      if(sqrt((i-f11_x)*(i-f11_x)+(j-f11_y)*(j-f11_y))+
				 sqrt((i-f21_x)*(i-f21_x)+(j-f21_y)*(j-f21_y))<2*a1/sky_size_arcsec*img_size)
				{
				  for(vector<double>::size_type k=0;k<freq_vec.size();++k)
				    {
				      double freq=freq_vec[k]*1E6;
				      img_vec[k](mod(i,img_size),mod(j,img_size))+=calc_Tb(fr1_lobe_spec(fr1.lobe_flux[0],freq),s_lobe1,freq);
				    }
				}
			    }
			}

		      xmin = (int)(x_lobe_2 - a2 / sky_size_arcsec * img_size);
		      xmax = (int)(x_lobe_2 + a2 / sky_size_arcsec * img_size);
		      ymin = (int)(y_lobe_2 - a2 / sky_size_arcsec * img_size);
		      ymax = (int)(y_lobe_2 + a2 / sky_size_arcsec * img_size);
		      for(int i = xmin; i <= xmax; ++i)
			{
			  for(int j = ymin; j <= ymax; ++j)
			    {
			      if(sqrt((i-f12_x)*(i-f12_x)+(j-f12_y)*(j-f12_y))+
				 sqrt((i-f22_x)*(i-f22_x)+(j-f22_y)*(j-f22_y))<2*a2/sky_size_arcsec*img_size)
				{
				  for(vector<double>::size_type k=0;k<freq_vec.size();++k)
				    {
				      double freq=freq_vec[k]*1E6;
				      img_vec[k](mod(i,img_size),mod(j,img_size))+=calc_Tb(fr1_lobe_spec(fr1.lobe_flux[1],freq),s_lobe2,freq);
				    }
				}
			    }
			}
		      
		      //
		      //if(sqrt((x_core-x_lobe_1)*(x_core-x_lobe_1)+(y_core-y_lobe_1)*(y_core-y_lobe_1))>=2)
		      {
			//   cout<<x_core<<"\t"<<y_core<<"\t"<<x_lobe_1<<"\t"<<y_lobe_1<<"\t"<<fr1.lobe_position_angle[0]<<endl;  
		      }
		      
		    }
		}
	    }
	  
	  old_galaxy=galaxy;
	}
    }
  
  for(size_t i=0; i<freq_vec.size(); ++i)
    {
      char *fname;
      cfitsfile ff;
      asprintf(&fname, "%s%.2f.fits", output_prefix, freq_vec[i]);
      ff.create(fname);
      ff << img_vec[i];
      printf("INFO: saved map '%s'\n", fname);
      free(fname);
    }
}
