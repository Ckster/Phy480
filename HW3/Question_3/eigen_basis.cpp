// The higher the dimension size, the better the convergence is. This is predictable
// since the real basis should be infinite. b=1.9 converged much more closely to the
//actual wavefunction than b=1.9 did. The higher b was made the flatter the curve got, but too small of a b
// (smaller than 1.9) made the curve too tall and skinny. There probably isnt a b that would
//make the plots converge, especially not with the basis sizes allowed. I would guess
// the plot that gives the least error is dimension 20 with a b between 1.9 and 2.5.



// Program to rebuild the wave function for a quantum harmonic oscillator
// experincing a coulomb potential using GSL eigenvaule/eigenvector routines.
// include files
#include <iostream>		// note that .h is omitted
#include <iomanip>		// note that .h is omitted
#include <cmath>
#include <fstream>
using namespace std;

#include <gsl/gsl_eigen.h>	        // gsl eigensystem routines
#include <gsl/gsl_integration.h>	// gsl integration routines

// structures and function prototypes
typedef struct			// structure holding Hij parameters
{
  int i;			// 1st matrix index
  int j;			// 2nd matrix index

  double mass;			// particle mass
  double b_ho;			// harmonic oscillator parameter
  int potential_index;		// indicates which potential to use
}
hij_parameters;

typedef struct			// structure holding potential parameters
{
  double param1;		// any three parameters
  double param2;
  double param3;
}
potential_parameters;

// potentials
double V_coulomb (double r, potential_parameters * potl_params_ptr);
double V_square_well (double r, potential_parameters * potl_params_ptr);
double V_morse (double r, potential_parameters * potl_params_ptr);

// i'th-j'th matrix element of Hamiltonian in ho basis
double Hij (hij_parameters ho_parameters);
double Hij_integrand (double x, void *params_ptr);

// harmonic oscillator routines from harmonic_oscillator.cpp
extern double ho_radial (int n, int l, double b_ho, double r);
extern double ho_eigenvalue (int n, int l, double b_ho, double mass);

// to square double precision numbers
inline double sqr (double x)  {return x*x;}

//************************** main program ***************************
int
main ()
{
  hij_parameters ho_parameters;  // parameters for the Hamiltonian


  ho_parameters.potential_index = 1; //Choose the coulomb potential




  double mass = 1;		 // measure mass in convenient units
  ho_parameters.mass = mass;






    ofstream out_1_19 ("eigen_basis_dim=1_b=1.9.dat");	// open the output file for dim 1 b=1.9
    ofstream out_5_19 ("eigen_basis_dim=5_b=1.9.dat"); // open the output file for dim 5 b=1.9
    ofstream out_10_19 ("eigen_basis_dim=10_b=1.9.dat"); // open the output file for dim 10 b=1.9
    ofstream out_20_19 ("eigen_basis_dim=20_b=1.9.dat"); // open the output file for dim 20 b=1.9
    ofstream out_1_25 ("eigen_basis_dim=1_b=2.5.dat");	// open the output file for dim 1 b=2.1
    ofstream out_5_25 ("eigen_basis_dim=5_b=2.5.dat"); // open the output file for dim 5 b=2.1
    ofstream out_10_25 ("eigen_basis_dim=10_b=2.5.dat"); // open the output file for dim 10 b=2.1
    ofstream out_20_25 ("eigen_basis_dim=20_b=2.5.dat"); // open the output file for dim 20 b=2.1
    int dimensions[]={1,5,10,20}; // These are the basis dimensions to be used in the sum
    double b_s[]={1.9,2.5}; // The two different osciallator paramaters that will be used
    double u = 0; //initialize the wavefunction at zero

    for (int k=0;k<2;k++){ //loop over the different b's
      double b_ho = b_s[k];
      ho_parameters.b_ho = b_ho;
    for (int i=0;i<4;i++){ //loop over the different dimensions
      int dimension = dimensions[i];

    for (double r = 0; r<15; r+=0.01){ //loop through values of r so the wavefunction can be plotted
      for (int j=0;j<dimension;j++){
        ho_parameters.i = 0; // we want the ground state wavefunction

    	  ho_parameters.j = j; // loop through different basis functions corresponding to each dimension
      u+=-Hij (ho_parameters)*ho_radial (1, 0, b_ho, r);  // n=0,l=0, for ground state. (n_0=1).
    }
    if (b_ho==1.9){ //output data to the correct data files. One for each b and 4 for each total basis dimension
    if (dimension==1){
    out_1_19 << scientific << setprecision (8)
      << r << "   "
			<< u/pow(b_ho,.5) << endl; // Normalize the wavefuntion! Divide by 1/sqrt(b).
    }
    if (dimension==5){
      out_5_19 << scientific << setprecision (8)
        << r << "   "
  			<< u/pow(b_ho,.5) << endl;
    }
    if (dimension==10){
      out_10_19 << scientific << setprecision (8)
        << r << "   "
  			<< u/pow(b_ho,.5) << endl;
    }
    if (dimension==20){
      out_20_19 << scientific << setprecision (8)
        << r << "   "
  			<<u/pow(b_ho,.5) << endl;
    }
  }
  if (b_ho==2.5){
  if (dimension==1){
  out_1_25 << scientific << setprecision (8)
    << r << "   "
    << u/pow(b_ho,.5) << endl;
  }
  if (dimension==5){
    out_5_25 << scientific << setprecision (8)
      << r << "   "
      << u/pow(b_ho,.5) << endl;
  }
  if (dimension==10){
    out_10_25 << scientific << setprecision (8)
      << r << "   "
      << u/pow(b_ho,.5) << endl;
  }
  if (dimension==20){
    out_20_25 << scientific << setprecision (8)
      << r << "   "
      << u/pow(b_ho,.5) << endl;
  }
}
      u = 0;

    }
  }
}



  out_1_19.close ();
  out_5_19.close ();
  out_10_19.close ();
  out_20_19.close ();
  out_1_25.close ();
  out_5_25.close ();
  out_10_25.close ();
  out_20_25.close ();
  return (0);			// successful completion
}

//************************************************************


//************************** Hij ***************************
//
// Calculate the i'th-j'th matrix element of the Hamiltonian
//  in a Harmonic oscillator basis.  This routine just passes
//  the integrand Hij_integrand to a GSL integration routine
//  (gsl_integration_qagiu) that integrates it over r from 0
//  to infinity
//
// Take l=0 only for now
//
//*************************************************************
double
Hij (hij_parameters ho_parameters)
{
 gsl_integration_workspace *work = gsl_integration_workspace_alloc (1000);
 gsl_function F_integrand;

 double lower_limit = 0.;	// start integral from 0 (to infinity)
 double abs_error = 1.0e-8;	// to avoid round-off problems
 double rel_error = 1.0e-8;	// the result will usually be much better
 double result = 0.;		// the result from the integration
 double error = 0.;		// the estimated error from the integration

 void *params_ptr;		// void pointer passed to function

 params_ptr = &ho_parameters;	// we'll pass i, j, mass, b_ho

 // set up the integrand
 F_integrand.function = &Hij_integrand;
 F_integrand.params = params_ptr;

 // carry out the integral over r from 0 to infinity
 gsl_integration_qagiu (&F_integrand, lower_limit,
			 abs_error, rel_error, 1000, work, &result, &error);
 // eventually we should do something with the error estimate

 return (result);		// send back the result of the integration
}

//************************** Hij_integrand ***************************
//
// The integrand for the i'th-j'th matrix element of the
//  Hamiltonian matrix.
//   * uses a harmonic oscillator basis
//   * the harmonic oscillator S-eqn was used to eliminate the
//      2nd derivative from the Hamiltonian in favor of the
//      HO energy and potential.  This was checked against an
//      explicit (but crude) 2nd derivative (now commented).
//
//************************************************************
double
Hij_integrand (double x, void *params_ptr)
{
 potential_parameters potl_params;	// parameters to pass to potential
 double Zesq;			// Ze^2 for Coulomb potential
 double R, V0;			// radius and depth of square well
 int potential_index;		// index 1,2,... for potental

 int l = 0;			// orbital angular momentum
 int n_i, n_j;			// principal quantum number (1,2,...)
 double mass, b_ho;		// local ho parameters
 double hbar = 1.;		// units with hbar = 1
 double omega;			// harmonic oscillator frequency
 double ho_pot;		// value of ho potential at current x

 // define variables for debugging 2nd derivative
 /*
 double h = 0.01;
 double fp, f, fm, deriv2;
 */

 n_i = ((hij_parameters *) params_ptr)->i + 1;	// n starts at 1
 n_j = ((hij_parameters *) params_ptr)->j + 1;
 mass = ((hij_parameters *) params_ptr)->mass;
 b_ho = ((hij_parameters *) params_ptr)->b_ho;
 omega = hbar / (mass * b_ho * b_ho);	// definition of omega
 ho_pot = (1. / 2.) * mass * (omega * omega) * (x * x);	// ho pot'l
 potential_index = ((hij_parameters *) params_ptr)->potential_index;

 // debugging code to calculate 2nd derivative by hand
 /*
 f = ho_radial (n_j, l, b_ho, x);
 fp = ho_radial (n_j, l, b_ho, x + h);
 fm = ho_radial (n_j, l, b_ho, x - h);
 deriv2 = -((fp - f) - (f - fm)) / (h * h) / (2. * mass);
 */

 // set up the potential according to potential index
 switch (potential_index)
   {
   case 1:			// coulomb
     Zesq = 1.;
     potl_params.param1 = Zesq;
     return (ho_radial (n_i, l, b_ho, x)
	      * (ho_eigenvalue (n_j, l, b_ho, mass) - ho_pot
		 + V_coulomb (x, &potl_params))
	      * ho_radial (n_j, l, b_ho, x));
     break;
   case 2:			// square well
     R = 1.;
     V0 = 50.;
     potl_params.param1 = V0;
     potl_params.param2 = R;
     return (ho_radial (n_i, l, b_ho, x)
	      * (ho_eigenvalue (n_j, l, b_ho, mass) - ho_pot
		 + V_square_well (x, &potl_params))
	      * ho_radial (n_j, l, b_ho, x));
     break;
   default:
     cout << "Shouldn't get here!\n";
     return (1);
     break;
   }


 // debugging code to use crude 2nd derivative
 // return (ho_radial (n_i, l, b_ho, x)
 //	  * (deriv2 + V_coulomb (x, &potl_params)
 //	     * ho_radial (n_j, l, b_ho, x)));

}

//************************** Potentials *************************

//************************** V_coulomb ***************************
//
// Coulomb potential with charge Z:  Ze^2/r
//  --> hydrogen-like atom
//
//   Zesq stands for Ze^2
//
//**************************************************************
double
V_coulomb (double r, potential_parameters * potl_params_ptr)
{
 double Zesq = potl_params_ptr->param1;

 return (-Zesq / r);
}

//**************************************************************

//************************* V_square_well **********************
//
// Square well potential of radius R and depth V0
//
//**************************************************************
double
V_square_well (double r, potential_parameters * potl_params_ptr)
{
 double V0 = potl_params_ptr->param1;
 double R = potl_params_ptr->param2;

 if (r < R)
   {
     return (-V0);		// inside the well of depth V0
   }
 else
   {
     return (0.);		// outside the well
   }
}

//************************************************************

//************************** V_morse ***************************
//
// Morse potential with equilibrium bond length r_eq and potential
//  energy for bond formation D_eq
//
//**************************************************************
double
V_morse (double r, potential_parameters * potl_params_ptr)
{
 double D_eq = potl_params_ptr->param1;
 double r_eq = potl_params_ptr->param2;

 return ( D_eq * sqr(1. - exp(-(r-r_eq))) );
}

//**************************************************************
