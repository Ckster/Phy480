extern double simpsons_rule ( int num_pts, double x_min, double x_max,
                       double (*integrand) (double x) ); //Simpsons Rule
extern double milne_rule ( int num_pts, double x_min, double x_max,
                      double (*integrand) (double x) );//Milne Rule
extern double my_gsl_rule(int i); //Gauss-Kronrod 21 point
//end: function prototypes
