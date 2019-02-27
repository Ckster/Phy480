///////////////////SLOPE-ANALYSIS////////////////////////////
// THE SLOPE OF THE LOG-LOG PLOT FOR SIMPSONS IS -3.99 WHEN FITTED FROM .5-1.3. THIS MEANS
//THAT IN THIS REGION THE ERROR SCALES LIKE ~ 1/N^4, AS IT SHOULD. BY LOOKING AT THE GRAPH
//THE MAX OPTIMAL NUMBER OF ITERATIONS SHOULD BE ~ 12.

//THE SLOPE OF THE LOG-LOG PLOT FOR MILNE IS -8.6 WHEN FITTED FROM 0-1.3. THIS MEANS
// THAT IN THIS REGION THE ERROR SCALES AS 1/N^8.6, AND IT SHOULD SCALE AS 1/N^8. HONING
// IN FURTHER ON THE CORRECT REGION COULD MAKE THIS CLOSER TO 8. BY LOOKING AT THE GRAPH THE
// OPTIMAL NUMBER OF ITERATIONS SHOULD BE ~ 4. AFTER THIS THE ERROR OSCILLATES RANDOMLY, BUT IS LOW.



#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "gsl/gsl_integration.h"// GSL integration routines for Gauss_Kronrod21.
using namespace std;

#include "integ_hw2_routines.h"	// integration routines for Gauss and Milne

float my_integrand (float x);
const double answer = 2.2955871493926;  //Answer to integral sqrt(1+x^2)

int
main ()
{
  // set up the integration specifiction
  const int max_intervals = 501;	// maximum number of intervals
  const double lower = -1.0;	// lower limit of integration
  const double upper = 1.0;	// upper limit of integration
  double result = 0.;  // approximate answer

  // open the output file stream
  ofstream integ_out_Simpsons ("integ_hw2_Simpsons.dat");	// save Simpson's data in integ_hw2_Simpsons.dat
  ofstream integ_out_milne ("integ_hw2_milne.dat"); //save Milne data into integ_hw2_milne.dat
  ofstream integ_out_GK21("integ_hw2_GK21.dat"); // save Gauss-Kronrod 21pt data in integ_hw2_GK21.dat
  integ_out_Simpsons << "#  N_Simpson  Simpsons  " << endl;
  integ_out_milne << "#  N_Milne  Milne " << endl;
  integ_out_GK21 << "#  N_GK21  GK21 " << endl;
  integ_out_Simpsons << "#-----------------------------------------" << endl;
  integ_out_milne << "#-----------------------------------------" << endl;
  integ_out_milne <<"#------------------------------------------" << endl;
  // Simpson's rule requires an odd number of intervals (2*i)+1 i=1,2,3....max_intervals
  for (int i = 3; i <= max_intervals; i += 2)
  {
    integ_out_Simpsons << setw(4) << log10(i);
    result = simpsons_rule (i, lower, upper, &my_integrand);
    integ_out_Simpsons << "  " << setprecision(10) << log10(fabs (((result - answer)/answer)));
    //output relative error of the integration for Simpsons Rule
    integ_out_Simpsons<<endl;
  }
 //Milne rule has odd number of intervals given by (4*i)+1. i =1,2,3..max_intervals
  for (int i = 5; i <= max_intervals; i += 4)
  {
    integ_out_milne << setw(4) << log10(i);


    result = milne_rule (i, lower, upper, &my_integrand);
    integ_out_milne << "  " << setprecision(10) << log10(fabs (((result - answer)/answer)));
    //output relative error of the integration for Milne Rule


    integ_out_milne << endl;
  }
   //Gauss_Kronrod 21 needs intervals spaced as (20*i)+1 i=1,2,3...max_intervals
    for (int i = 21; i <= max_intervals; i += 20)
    {
    integ_out_GK21 << setw(4) << log10(i);
    result = my_gsl_rule(i);
    integ_out_GK21 << "  " << setprecision(10) << log10(fabs (((result - answer)/answer)));
    //output relative error of the integration for Gauss_Kronrod 21 point
    integ_out_GK21 << endl;
  }





  cout << "Simpson's data stored in integ_hw2_Simpsons.dat\n";
  cout << "Milne data stored in integ_hw2_Milne.dat\n";
  cout << "Gauss-Kronrod 21 data stored in integ_hw2_GK21.dat\n";
  integ_out_Simpsons.close ();
  integ_out_milne.close ();
  integ_out_GK21.close ();
  return (0);
}

//************************************************************************

// the function we want to integrate
float
my_integrand (float x)
{
  return (sqrt(1+(x*x)));
}
