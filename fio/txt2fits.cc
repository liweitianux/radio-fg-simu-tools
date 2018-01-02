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
  //To see how many arguments the user supplies
  if(argc!=3)
    {
      cerr<<argv[0]<<" <input text file> <output fits file>"<<endl;
      return -1;
    }

  //define an object that can read data from a text file from disk
  ifstream ifs(argv[1]);

  //define the varible to store the size of the image
  int height,width;
  //read the size from the text file
  ifs>>height>>width;
  
  //use the size to initial a matrix, its name is img
  Array<double,2> img(height,width);
  
  //go on to read the following data from the file

  for(int i=0;i!=img.extent(0);++i)
    {
      for(int j=0;j!=img.extent(1);++j)
	{
	  ifs>>img(i,j);
	}
    }

  //define a fitsfile object
  cfitsfile ff;
  //create the fitsfile
  ff.create(argv[2]);
  //trans the data from the matrix to the fits file.
  ff<<img;
  //done
}
