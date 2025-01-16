// ===========================================================================
// FUNCTION FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim, int degree>
void customPDE<dim,degree>::setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC){
    // ---------------------------------------------------------------------
    // ENTER THE INITIAL CONDITIONS HERE
    // ---------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the initial condition for each variable
    // according to its variable index
      double center[3] = {10,12,0.0};
      double radius = userInputs.domain_size[0]/10.0;

      double dist;
      dist = 0.0;
      for (unsigned int dir = 0; dir < dim; dir++){
          dist += (p[dir]-center[dir])*(p[dir]-center[dir]);
      }
      dist = std::sqrt(dist); 
 
     scalar_IC = 0;
        if (index == 0){
                if((std::abs(p[1]-userInputs.domain_size[1]/2.0 + dx*0.5) < dx) && (p[0] < clength)) {
                        scalar_IC = 1.0;
                }

        }

      if (index ==1){
          for (unsigned int d=0; d<dim; d++){
              vector_IC(d) = 0.0;
          }
      }
      if (index == 2) scalar_IC = 0.0;
      if (index == 3) scalar_IC = 1.0;
      if (index == 4) scalar_IC = 1.0;
      // --------------------------------------------------------------------------
  }

  // ===========================================================================
  // FUNCTION FOR NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS
  // ===========================================================================

  template <int dim, int degree>
  void customPDE<dim,degree>::setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC)
  {
      // --------------------------------------------------------------------------
      // ENTER THE NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS HERE
      // --------------------------------------------------------------------------
      // Enter the function describing conditions for the fields at point "p".
      // Use "if" statements to set the boundary condition for each variable
      // according to its variable index. This function can be left blank if there
      // are no non-uniform Dirichlet boundary conditions. For BCs that change in
      // time, you can access the current time through the variable "time". The
      // boundary index can be accessed via the variable "direction", which starts
      // at zero and uses the same order as the BC specification in parameters.in
      // (i.e. left = 0, right = 1, bottom = 2, top = 3, front = 4, back = 5).
      // -------------------------------------------------------------------------
    if (index == 1){
	double x = (p[0] - vel_nom*time - clength);
	double y = p[1] - userInputs.domain_size[1]/2.0+dx*0.5;
	double r = std::sqrt(x*x+y*y);
	double pi = 3.14159265359;
	double theta = std::atan2(y,x);
        double mu = CIJ_base[dim][dim];
	double lambda = CIJ_base[0][0] - 2.0*mu;
	double nu = lambda/2.0/(lambda+mu);
	double kappa = 3.0-4.0*nu;
    	vector_BC[0] = 0.5*(KI_nom/mu)*std::sqrt(0.5*r/pi)*std::cos(0.5*theta)*(kappa - std::cos(theta));
    	vector_BC[1] = 0.5*(KI_nom/mu)*std::sqrt(0.5*r/pi)*std::sin(0.5*theta)*(kappa - std::cos(theta));
    }
  }
