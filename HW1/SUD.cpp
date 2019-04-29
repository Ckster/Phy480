//Compile using g++ -o SUD.x SUC.cpp
//The sum down method is more accurate because smaller numbers are being
//added to larger numbers, and in the forward method large numbers are being
//added to small numbers and roundoff errors occur, making it less accurate.
//When you are adding two numbers, then, you want to add from the smallest
//number to the biggest number. For large values of N the errors begin to oscillate
// and for small values of N it is a power law.



#include <iostream>	     // this has the cout, cin definitions
#include <iomanip>
using namespace std;     // if omitted, then need std::cout, std::cin
#include <fstream>
#include <cmath>

int main()
{
  int N;

  float diff;
  float sum_up = 0;  // set the initial sum to zero
  float sum_down = 0;
  float n=1;
  
  //INPUT FOR N SHOULD BE LARGE ~10^7


  ofstream sum_out ("sum_out.dat");	// save data in sum_out.dat
  sum_out << "#  N   Sum Up  Sum Down Difference " << endl;
  sum_out << "#----------------------------------" << endl;
  while(n<100000){
      N=n;
      sum_up = 0;
      for (float i = 1; i <= N; i+=1){

        sum_up  += 1.0/i;

          }
      sum_down = 0;
      for (int i = N; i >= 1; i -= 1){

        sum_down += 1./i;
        }
      diff = (abs(sum_up - sum_down)/((1./2.) * (abs(sum_up) + abs(sum_down))));
      sum_out << setw(4) << n;
      sum_out << "  " << setprecision(12) << sum_up;
      sum_out << "  " << setprecision(12) << sum_down;
      sum_out << "  " << setprecision(12) << diff;
      sum_out << endl;
      n+=1;
  }
cout << "data stored in sum.dat\n";
sum_out.close ();
return 0;

}
