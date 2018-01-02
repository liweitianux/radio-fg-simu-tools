#include <iostream>
#include <fstream>
#include <random/normal.h>
using namespace std;
using namespace ranlib;

int n=0;
int main(int argc,char* argv[])
{
  if (argc!=2)
    {
      cerr<<argv[0]<<"<cluster number> "<<endl;
      return 0;
    }
  int N=(int)atof(argv[1]);

  Normal<double>(0,1).seed((unsigned int)time(0));
   
  double beta0=0.66;
  for (int i=0;;++i)
    {
    repeat:
     
      Normal<double> rnd(beta0,0.33);
      double beta=rnd.random();
      if (beta>=0.33&&beta<=0.99&&n<N)
	{
	 
	  cout<<beta<<endl;
	  n++;	
	}
      else 
	{
	  if (n<N)
	    {
	      goto repeat;
	    }
	  else return 0;
	}
    }
 
}
