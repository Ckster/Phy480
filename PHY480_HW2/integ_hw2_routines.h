extern double simpsons_rule ( int num_pts, float x_min, float x_max,
                       float (*integrand) (float x) ); //Simpsons Rule
extern double milne_rule ( int num_pts, float x_min, float x_max,
                      float (*integrand) (float x) );//Milne Rule
extern double my_gsl_rule(int i); //Gauss-Kronrod 21 point
//end: function prototypes 
