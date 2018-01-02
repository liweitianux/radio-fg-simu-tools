#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

const double k=-2.1568;
const double b=3.45736;
const double delta_b=0.8;

int main(int argc,char* argv[])
{
  for (int i=-10;i<20;++i)
    {
      double x=pow(10,i/10.);
      double y=k*log10(x)+b;
      cout<<x<<"\t0\t"<<y<<"\t0"<<endl; 
    }
  cout<<"no no no"<<endl;
  for (int i=-10;i<20;++i)
    {
      double x=pow(10,i/10.);
      double y=k*log10(x)+b-delta_b;
      cout<<x<<"\t0\t"<<y<<"\t0"<<endl; 
    }
  cout<<"no no no"<<endl;
  for (int i=-10;i<20;++i)
    {
      double x=pow(10,i/10.);
      double y=k*log10(x)+b+delta_b;
      cout<<x<<"\t0\t"<<y<<"\t0"<<endl; 
    }
  cout<<"no no no"<<endl;
  cout<<"3  0  1  0"<<endl;
  cout<<"20  0  1  0"<<endl;
}

