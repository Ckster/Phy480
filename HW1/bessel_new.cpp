//When the error is close to one it means the difference of values in up
//and down is large. When X is near 10 they are both seen to be accurate, the
//error is small. From x=1 to x=10, there is a power law and for x>50 the upward
//method is better. 



// include files
#include <iostream>		// note that .h is omitted
#include <iomanip>		// note that .h is omitted
#include <fstream>		// note that .h is omitted
#include <cmath>
#include <gsl/gsl_sf_bessel.h>
using namespace std;		// we need this when .h is omitted

// function prototypes
double down_recursion (double x, int n, int m);	// downward algorithm
double up_recursion (double x, int n);	        // upward algorithm


// global constants
const double xmax = 100.0;	// max of x
const double xmin = 0.1;	// min of x >0
const double step = 0.1;	// delta x
const int order = 10;		// order of Bessel function
const int start = 50;		// used for downward algorithm

//********************************************************************
int
main ()
{
  double ans_down, ans_up, relative_diff, gsl;

  // open an output file stream
  ofstream my_out ("bessel_new.dat");

  my_out << "# Spherical Bessel functions via up and down recursion   Relative Difference"
         << endl;

  // step through different x values
  for (double x = xmin; x <= xmax; x += step)
    {
      ans_down = down_recursion (x, order, start);
      ans_up = up_recursion (x, order);
      relative_diff = abs(ans_down - ans_up)/(abs(ans_up) + abs(ans_down));
      gsl = gsl_sf_bessel_J1 (x);
      my_out << " " << setprecision (6) << x;
	    my_out << setw (13) << ans_down;
    	my_out << setw (13) << ans_up;
      my_out << setw (13) << gsl;
      my_out << " " << setprecision (6) << relative_diff
        << endl;
    }
  cout << "data stored in bessel_new.dat." << endl;

  // close the output file
  my_out.close ();
  return (0);
}


//------------------------end of main program-----------------------

// function using downward recursion
double
down_recursion (double x, int n, int m)
{
  double j[start + 2];		// array to store Bessel functions
  j[m + 1] = j[m] = 1.;		// start with "something" (choose 1 here)
  for (int k = m; k > 0; k--)
    {
      j[k - 1] = ((2.* double(k) + 1.) / x) * j[k] - j[k + 1];  // recur. rel.
    }
  double scale = (sin (x) / x) / j[0];	// scale the result
  return (j[n] * scale);
}


//------------------------------------------------------------------

// function using upward recursion
double
up_recursion (double x, int n)
{
  double term_three = 0.;
  double term_one = (sin (x)) / x;	// start with lowest order
  double term_two = (sin (x) - x * cos (x)) / (x * x);	// next order
  for (int k = 1; k < n; k += 1)	// loop for order of function
    { // recurrence relation
      term_three = ((2.*double(k) + 1.) / x) * term_two - term_one;
      term_one = term_two;
      term_two = term_three;
    }
  return (term_three);
}
