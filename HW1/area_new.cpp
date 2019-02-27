//  file: area.cpp
//
//  This program calculates the area of a circle, given the radius.
//
//  Programmer:  Dick Furnstahl  furnstahl.1@osu.edu
//
//  Revision history:
//      02-Jan-2004  original version, for 780.20 Computational Physics
//      01-Jan-2010  updates to "To do" wishlist
//      12-Jan-2016  comment out "using namespace std;"
//
//  Notes:
//   * compile with:  "g++ -o area.x area.cpp"
//
//  To do:
//   1. output the answer with higher precision (more digits)
//   2. use a "predefined" value of pi or generate it with atan
//   3. define an inline square function
//   4. split the calculation off into a function (subroutine)
//   5. output to a file (and/or input from a file)
//   6. add checks of the input (e.g., for non-positive radii)
//   7. rewrite using a Circle class
//
//*********************************************************************//

// include files
#include <iostream>	     // this has the cout, cin definitions
#include<iomanip>        // this includes the "setprecision"
#include<math.h>         // include this to get Pi as a built-in package
#include<fstream>
using namespace std;     // if omitted, then need std::cout, std::cin

//*********************************************************************//

//const double pi = 3.1415926535897932385;   // define pi as a constant
inline double square(double a)   //defined inline function to square radius and return variable
{
  return (a * a);
}
inline double areacalc(double z) // defined inline function to calculate area using output from square
{
  return (M_PI * z);
}
int
main ()
{

  double radius;    // every variable is declared as int or double or ...
  double z;        // have to declare the output of square; z, as double
  ofstream area_out ("area.dat");
  area_out << " Radius  Area" << endl;
  cout << "Enter the radius of a circle: ";	// ask for radius
  cin >> radius;
  z = square(radius); // define z as the output of square(radius)
  //double area = M_PI * square(radius);	// standard area formula
  if (radius < 0){
    cout<< "Radius must be greater than zero";
    return 0;

  }
  cout << "radius = " << radius << ",  area = "<<setprecision(12)<< areacalc(z);
  cout << "data stored in area.dat\n";
  area_out << "  " << setprecision(4) << radius;
  area_out << "  " << setprecision(12) << areacalc(z);
  area_out.close ();

  return 0;			// "0" for successful completion
}

//*********************************************************************//
