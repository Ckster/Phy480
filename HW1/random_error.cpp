#include <iostream>	     // this has the cout, cin definitions
#include <iomanip>
using namespace std;     // if omitted, then need std::cout, std::cin
#include <fstream>
#include<math.h>

int
main ()
{
  ofstream my_file; // add this to initialize output data file
  my_file.open("random.dat");
  float x_f = 0. ;
  double x_d = 0.;
  double error;

  for (int i = 1; i < 10000; i+=1){  //create a loop to take the square root of a bunch of numbers
    x_f=sqrt(float(i));
    x_d=sqrt(double(i));
    error = ((x_d - x_f)/(x_f));   // take the relative error but keep the sign, so on abs
    my_file<<" "<<i<<" " << x_d << " " << x_f << " " <<setprecision(16)<< error << " " << endl;
  }
  my_file.close();
  return 0;



}
