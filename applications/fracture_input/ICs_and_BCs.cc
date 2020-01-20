// ===========================================================================
// FUNCTION FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim, int degree>
void customPDE<dim,degree>::setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC, const std::vector<double> &data){
    // ---------------------------------------------------------------------
    // ENTER THE INITIAL CONDITIONS HERE
    // ---------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the initial condition for each variable
    // according to its variable index

	  double center[1][3] = {{0.5,0.5,0.5}};
	  double rad[1] = {userInputs.domain_size[0]/40.0};
	  double dx=userInputs.domain_size[0]/((double) userInputs.subdivisions[0])/std::pow(2.0,userInputs.refine_factor);
	  double dist;
	  
	  scalar_IC = 0;

	  //for (unsigned int i=0; i<1; i++){
		//  dist = 0.0;
		//  for (unsigned int dir = 0; dir < dim; dir++){
		//	  dist += (p[dir]-center[i][dir]*userInputs.domain_size[dir])*(p[dir]-center[i][dir]*userInputs.domain_size[dir]);
		//  }
		//  dist = std::sqrt(dist);

      	if (index == 0){
		if((std::abs(p[1]-userInputs.domain_size[1]/2.0) < cwidth*0.5)&&(p[0] < clength)) {
			scalar_IC = 1.0;
		}
		
	}

	//}

      if (index ==1){
          for (unsigned int d=0; d<dim; d++){
              vector_IC(d) = 0.0;
          }
      }
      if (index > 2){
         int i = (int)std::floor(p[0]*(128.0)/userInputs.domain_size[0]);
         int j = (int)std::floor(p[1]*(128.0)/userInputs.domain_size[1]);
	scalar_IC = data[j*257+i];
        if (index == 4) scalar_IC = 1.0;
      }

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
        if (direction == 3){
            vector_BC[0]=0.0;
	    if(stepped_strain == false){
	 // working simple BC
                vector_BC[1]=u_init+u_step*time/step_t;
	    }else{
         // BC to isolate viscous effects
	        vector_BC[1]=u_init+u_step*std::ceil((time - userInputs.dtValue)/step_t);
	    }
        }
    }

  }
