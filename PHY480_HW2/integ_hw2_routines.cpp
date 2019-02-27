#include <cmath>
#include "integ_hw2_routines.h"   // integration routine simspson and Milne
#include "gsl/gsl_integration.h" // gsl integration routines
double my_gsl_integrand (double x, void *params_ptr);
double simpsons_rule ( int num_pts, float x_min, float x_max,
                      float (*integrand) (float x) )


{
   double interval = ((x_max - x_min)/float(num_pts - 1));  // called h in notes
   double sum=  0.;  // initialize integration sum to zero

   for (int n=2; n<num_pts; n+=2)                // loop for odd points
   {
     double x = x_min + interval * float(n-1);
     sum += (4./3.)*interval * integrand(x);
   }
   for (int n=3; n<num_pts; n+=2)                // loop for even points
   {
      double x = x_min + interval * float(n-1);
     sum += (2./3.)*interval * integrand(x);
   }
   // add in the endpoint contributions
   sum +=  (interval/3.) * (integrand(x_min) + integrand(x_max));

   return (sum);
}
double milne_rule (int num_pts, float x_min, float x_max,
                      float (*integrand) (float x))
{

  double interval = ((x_max - x_min)/float(num_pts - 1));  // called h in notes
  double sum=  0.;  // initialize integration sum to zero

  for (int n=2; n<num_pts; n+=2)                // loop for odd points
  {
    double x = x_min + interval * float(n-1);
    sum += (64./45.)*interval * integrand(x);
  }
  for (int n=3; n<num_pts; n+=4)                // loop for first set of even points
  {
     double x = x_min + interval * float(n-1);
    sum += (24./45.)*interval * integrand(x);
  }
  for (int n=5; n<num_pts; n+=4)                // loop for second set of even points
  {
     double x = x_min + interval * float(n-1);
    sum += (28./45.)*interval * integrand(x);
  }
  // add in the endpoint contributions
  sum +=  (14./45.)*interval * (integrand(x_min) + integrand(x_max));

  return (sum);
}
double my_gsl_rule(int i) //initialize gsl integration function

  {

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (i); //creates workspace of
    //to hold i double precision intervals, integration results and their errors


    double result1, error; //declare result and error for gsl routine

    double lower_limit = 0.; //lower limit of integration

    double upper_limit = 1.; //upper limit of integration

    double alpha = 1.; //paramater for integrand

    void *params_ptr; //pointer to hold all paramaters (just alpha in this case)

    params_ptr = &alpha;

    gsl_function F; //input function for the gsl routine

    F.params = params_ptr; //assigns paramaters to F

    F.function = &my_gsl_integrand; //defines the integrand for F




    gsl_integration_qags (&F, lower_limit, upper_limit, 0, 1e-2, i,w, &result1, &error);
    //qags applies Gauss-Kronrod 21 point adaptively to get numerical result for integral
    gsl_integration_workspace_free (w); //frees the memory assosciated with workspace w.

    return(result1); //return nuerical result


  }




  double

 my_gsl_integrand (double x, void *params_ptr)  // create function that returns integrand for gsl

 {

   double alpha;

   alpha = *(double *) params_ptr; //Calls in paramater alpha


  return (alpha*sqrt(1+(x*x))); //Returns actual integrand
  }
