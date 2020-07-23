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
        //    dist += (p[dir]-center[i][dir]*userInputs.domain_size[dir])*(p[dir]-center[i][dir]*userInputs.domain_size[dir]);
        //  }
        //  dist = std::sqrt(dist);

   //     if (index == 0){
    //    if((std::abs(p[1]-userInputs.domain_size[1]/2.0) < cwidth*0.5)&&(p[0] < clength)) {
    //        scalar_IC = 1.0;
    //    }
     /*   if (index == 0){
		double y = p[1]-userInputs.domain_size[1]/2.0;
		double dist = std::sqrt(p[0]*p[0]+y*y);
		if (dist < std::sqrt(2.0*ell2)){ 
			scalar_IC = (1.0 - dist/std::sqrt(2.0*ell2))*(1.0 - dist/std::sqrt(2.0*ell2));
		}
        }*/
        

      if (index ==1){
          for (unsigned int d=0; d<dim; d++){
              vector_IC(d) = 0.0;
          }
      }

      // interpolation from regularly gridded input data
      if (index > 20){
          int r1[dim];
          int r2[dim];
          double dist[dim];
          for (int j=0;j<dim;++j){
          	  double dndx = (userInputs.bindata_size[j]-1.0)/userInputs.domain_size[j];
        	  r1[j] = (int) std::floor(p[j]*dndx);
        	  double rem = p[j] - r1[j]/dndx;
        	  r2[j] = (int) std::min(r1[j]+1,userInputs.bindata_size[j]-1);
        	  dist[j] = rem*dndx;
          }
          
          if (dim == 2){
             int ny = userInputs.bindata_size[1];
             double xint1 = data[r1[0]*ny+r1[1]]*(1.0-dist[0])  +
             			    data[r2[0]*ny+r1[1]]*dist[0]; 
             double xint2 = data[r1[0]*ny+r2[1]]*(1.0-dist[0]) +
             			    data[r2[0]*ny+r2[1]]*dist[0]; 
             scalar_IC = xint1*(1.0-dist[1]) + xint2*dist[1];
          
    	 }
   		 if (dim == 3){
   		     int ny = userInputs.bindata_size[1];
   		     int nz = userInputs.bindata_size[2];
	  //   std::cout << nz << " " << nz*ny*nz << " " << r1[0]*ny*nz + r1[1]*ny + r1[2] << "\n";
             double xint00 = data[r1[0]*ny*nz + r1[1]*ny + r1[2]]*(1.0-dist[0])  +
             			     data[r2[0]*ny*nz + r1[1]*ny + r1[2]]*dist[0]; 
	  //   std::cout << data[r1[0]*ny*nz + r1[1]*ny + r1[2]] << " " << data[r2[0]*ny*nz + r1[1]*ny + r1[2]] << " " << xint00 << "\n";
             double xint01 = data[r1[0]*ny*nz + r1[1]*ny + r2[2]]*(1.0-dist[0]) +
             			     data[r2[0]*ny*nz + r1[1]*ny + r2[2]]*dist[0] ; 
             double xint10 = data[r1[0]*ny*nz + r2[1]*ny + r1[2]]*(1.0-dist[0]) +
             			     data[r2[0]*ny*nz + r2[1]*ny + r1[2]]*dist[0] ; 
             double xint11 = data[r1[0]*ny*nz + r2[1]*ny + r2[2]]*(1.0-dist[0]) +
             			     data[r2[0]*ny*nz + r2[1]*ny + r2[2]]*dist[0] ;
             double yint0 = xint00*(1.0-dist[1]) + xint10*dist[1];
             double yint1 = xint01*(1.0-dist[1]) + xint11*dist[1];
             scalar_IC = yint0*(1.0-dist[2]) + yint1*dist[2];
	//     std::cout << r1[0] << " " << r1[1] << " " << r1[2] << " " << r2[0] << " " << r2[1] << " " << r2[2] << " "
	//	 << xint00 << " " << xint01 << " " << xint10 << " " << xint11 << " " << yint0 << " " << yint1 << " " << scalar_IC << "\n"; 		 
    	}
    }

        if (index == 3) scalar_IC = 1.0;
     /*   if (index == 3) {
		double distx = p[0]-userInputs.domain_size[0]/2.0;
		double disty = p[1]-userInputs.domain_size[1]/2.0;
		double rad = std::sqrt(distx*distx+disty*disty);
		scalar_IC = 0.5*(1.0+std::tanh(3.0*(rad-1.0)));
		if (scalar_IC < 1e-3) scalar_IC = 1e-3;
	}*/
       if (index == 4) scalar_IC = 1.0;
        if (index == 0) {
		scalar_IC = 0.0;
		double distx = p[0]-userInputs.domain_size[0]/2.0;
		double disty = p[1]-userInputs.domain_size[1]/2.0;
		double rad = std::sqrt(distx*distx+disty*disty);
		double pi = std::atan(1.0)*4.0;
		if (rad < 1.0) scalar_IC = 0.5*(1.0+std::cos(pi*rad));
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
/*        if (direction == 3){
            vector_BC[0]=0.0;
        if(stepped_strain == false){
     // working simple BC
                vector_BC[1]=u_init+u_step*time/step_t;
        }else{
         // BC to isolate viscous effects
            vector_BC[1]=u_init+u_step*std::ceil((time - userInputs.dtValue)/step_t);
        }
        }
*/
	// surfing BCs inspired by Brach et al. 2019
    }


  }
