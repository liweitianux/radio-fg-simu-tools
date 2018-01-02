#include <iostream>
#include <cmath>
#include <cstdlib>
#include <utility>
#include <fstream>
#include <cassert>
#include <vector>
#include "fio.h"

using namespace std;
using namespace ::blitz;



int main(int argc,char* argv[])
{
  //To see how many arguments the user supplies, 
  
  if(argc!=3)
    {
      cerr<<argv[0]<<" <input fits file> <output text file>"<<endl;
      return -1;
    }
  //define an cfitsfile object, its name is ff
  cfitsfile ff;
  //open the fits file with the name corresponding to 
  //the 1st argument
  ff.open(argv[1]);
  
  //define an 2-D matrix, its name is img
  Array<double,2> img;
  //transfer the data from the fitsfile to img
  ff>>img;
  
  //define an object that can be used to write text file to disk
  ofstream ofs(argv[2]);

  //firstly output the size of the image
  ofs<<img.extent(0)<<"\t"<<img.extent(1)<<endl;;
  
  //iterator over all the elements in the matrix, and output to the file
  for(int i=0;i!=img.extent(0);++i)
    {
      for(int j=0;j!=img.extent(1);++j)
	{
	  ofs<<img(i,j)<<" ";
	}
      ofs<<endl;
    }
  //done
}
