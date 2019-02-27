#include <iostream>	     // this has the cout, cin definitions
#include <iomanip>
using namespace std;     // if omitted, then need std::cout, std::cin
#include <fstream>

int main()
{
  int N;
  int i;
  float diff;
  float sum_up = 0;  // set the initial sum to zero
  float sum_down = 0;
  cout << "Enter N; the number of iterations: ";	// ask for N
  cin >> N;
  ofstream sum_out ("sum_out.dat");	// save data in sum_out.dat
  sum_out << "#  N   Sum Up  Sum Down Difference " << endl;
  sum_out << "#----------------------------------" << endl;
  for (int n = 1; n<=N; n++){
      sum_up = 0;
      for (int i = 1; i <= n; ++i){

        sum_up  += 1./i;

          }
      sum_down = 0;
      for (int i = n; i >= 1; i -= 1){

        sum_down += 1./i;
        }
      diff = (abs(sum_up - sum_down)/((1./2.) * (abs(sum_up) + abs(sum_down))));
      sum_out << setw(4) << n;
      sum_out << "  " << setprecision(12) << sum_up;
      sum_out << "  " << setprecision(12) << sum_down;
      sum_out << "  " << setprecision(12) << diff;
      sum_out << endl;
  }
cout << "data stored in sum.dat\n";
sum_out.close ();
return 0;

}
