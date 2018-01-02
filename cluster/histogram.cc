#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
double beta;
double data1,data2,data3,data4,data5,data6;
using namespace std;
int main(int argc,char* argv[])
{
  if(argc!=2)
    {
      cerr<<argv[0]<<" <data> "<<endl;
      return 0;
    }
  int total_n=0;
  double total_beta=0;
  ifstream ifs(argv[1]);
  const double resolution=0.05; 
  const double beta_max=5;
  vector<int> n(int(beta_max/resolution));
  while (1)
    {
       ifs>>beta;
      if (!ifs.good())
	{
	  break;
	}

      //    ifs>>data2>>data3>>data4>>data5>>data6>>beta;
     total_beta+=beta;
      total_n+=1;
      for (int i=0;i<beta_max/resolution;++i)
	{
	  if (beta>=i*resolution&&beta<(i+1)*resolution)
	    {
	      n[i]+=1;
	    }
	}
    }
      cout<<"!ave_beta="<<total_beta/total_n<<endl;
  for (int i=0;i<beta_max/resolution;++i)
    {
      cout<<i*resolution+resolution/2<<"\t"<<n[i]<<endl;
    }
}
