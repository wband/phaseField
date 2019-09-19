// List of variables and residual equations for the mechanics example application

// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadVariableAttributes(){

// Variable 0 - phase field
	set_variable_name				(0,"n");
	set_variable_type				(0,SCALAR);
	set_variable_equation_type		(0,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(0, "n, dndt");
    set_dependencies_gradient_term_RHS(0, "");
	
// Variable 1 - displacement field
	set_variable_name				(1,"u");
	set_variable_type				(1,VECTOR);
	set_variable_equation_type		(1,TIME_INDEPENDENT);

    set_dependencies_value_term_RHS(1, "");
    set_dependencies_gradient_term_RHS(1, "n, grad(u)");
    set_dependencies_value_term_LHS(1, "");
    set_dependencies_gradient_term_LHS(1, "n, grad(change(u))");

// Variable 2 - phase field update
	set_variable_name				(2,"dndt");
	set_variable_type				(2,SCALAR);
	set_variable_equation_type		(2,AUXILIARY);

    set_dependencies_value_term_RHS(2, "n, grad(u)");
    set_dependencies_gradient_term_RHS(2, "grad(n)");
	
}

// =============================================================================================
// explicitEquationRHS (needed only if one or more equation is explict time dependent)
// =============================================================================================
// This function calculates the right-hand-side of the explicit time-dependent
// equations for each variable. It takes "variable_list" as an input, which is a list
// of the value and derivatives of each of the variables at a specific quadrature
// point. The (x,y,z) location of that quadrature point is given by "q_point_loc".
// The function outputs two terms to variable_list -- one proportional to the test
// function and one proportional to the gradient of the test function. The index for
// each variable in this list corresponds to the index given at the top of this file.

template <int dim, int degree>
void customPDE<dim,degree>::explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {
// The phase field and the variable containing its update
scalarvalueType n = variable_list.get_scalar_value(0);
scalarvalueType dndt = variable_list.get_scalar_value(2);


// enforcing irreversibility: phase field must increase monotonically
for (unsigned int j=0; j<dndt.n_array_elements; ++j){
    if (dndt[j] > 0.0) dndt[j] = 0.0;
}

// update to phase field 
scalarvalueType eq_n = (n-constV(userInputs.dtValue*MnV)*dndt);

variable_list.set_scalar_value_term_RHS(0,eq_n);
}


// =============================================================================================
// nonExplicitEquationRHS (needed only if one or more equation is time independent or auxiliary)
// =============================================================================================
// This function calculates the right-hand-side of all of the equations that are not
// explicit time-dependent equations. It takes "variable_list" as an input, which is
// a list of the value and derivatives of each of the variables at a specific
// quadrature point. The (x,y,z) location of that quadrature point is given by
// "q_point_loc". The function outputs two terms to variable_list -- one proportional
// to the test function and one proportional to the gradient of the test function. The
// index for each variable in this list corresponds to the index given at the top of
// this file.

template <int dim, int degree>
void customPDE<dim,degree>::nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

 // --- Getting the values and derivatives of the model variables ---

//n
scalarvalueType n = variable_list.get_scalar_value(0);
scalargradType nx = variable_list.get_scalar_gradient(0);
//u
vectorgradType ux = variable_list.get_vector_gradient(1);


//compute strain tensor
dealii::VectorizedArray<double> E[dim][dim], S[dim][dim];
 for (unsigned int i=0; i<dim; i++){
 for (unsigned int j=0; j<dim; j++){
          E[i][j]= constV(0.5)*(ux[i][j]+ux[j][i]);
 }
 }

//compute stress tensor
computeStress<dim>(CIJ_Mg, E, S);


//compute RHS for mechanics solver (multiply stress by interpolation function)
 vectorgradType eqx_u;
for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
		eqx_u[i][j] = -S[i][j]*((constV(1.0)-2.0*n+n*n)+constV(k_small));
	}
}

// compute elastic energy density (without interpolation function)
dealii::VectorizedArray<double> elastic_energy=constV(0.0);

for (unsigned int i=0; i<dim; i++){
          for (unsigned int j=0; j<dim; j++){
                  elastic_energy += constV(0.5)*S[i][j]*E[i][j];
          }
}

// phase field update: (interpolation f'n)*(elastic energy)+bulk energy
scalarvalueType eq_n = (constV(2.0)*(n-constV(1.0))*elastic_energy + gc_ell*n);
// phase field update: gradient term (Laplacian in strong form)
scalargradType eqx_n = (constV(gc_ell*ell2)*nx);

// --- Submitting the terms for the governing equations ---
// mechanics
variable_list.set_vector_gradient_term_RHS(1,eqx_u);
// phase field
variable_list.set_scalar_value_term_RHS(2,eq_n);
variable_list.set_scalar_gradient_term_RHS(2,eqx_n);

}

/// =============================================================================================
// equationLHS (needed only if at least one equation is time independent)
// =============================================================================================
// This function calculates the left-hand-side of time-independent equations. It
// takes "variable_list" as an input, which is a list of the value and derivatives of
// each of the variables at a specific quadrature point. The (x,y,z) location of that
// quadrature point is given by "q_point_loc". The function outputs two terms to
// variable_list -- one proportional to the test function and one proportional to the
// gradient of the test function -- for the left-hand-side of the equation. The index
// for each variable in this list corresponds to the index given at the top of this
// file. If there are multiple elliptic equations, conditional statements should be
// sed to ensure that the correct residual is being submitted. The index of the field
// being solved can be accessed by "this->currentFieldIndex".

template <int dim, int degree>
void customPDE<dim,degree>::equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

// --- Getting the values and derivatives of the model variables ---
// phase field
scalarvalueType n = variable_list.get_scalar_value(0);
// displacement
vectorgradType Dux = variable_list.get_change_in_vector_gradient(1);

// --- Setting the expressions for the terms in the governing equations ---

vectorgradType eqx_Du;


// compute stiffness tensor
dealii::Tensor<2, CIJ_tensor_size, dealii::VectorizedArray<double> > CIJ;

if (n_dependent_stiffness == true){
for (unsigned int i=0; i<2*dim-1+dim/3; i++){
          for (unsigned int j=0; j<2*dim-1+dim/3; j++){
                  CIJ[i][j] = CIJ_Mg[i][j]*((constV(1.0)-2.0*n+n*n)+constV(k_small));
          }
}
}

//compute strain tensor
dealii::Tensor<2, dim, dealii::VectorizedArray<double> > E;

// couldn't use symmetrize before because of type issues
E = symmetrize(Dux);

//compute stress tensor
computeStress<dim>(CIJ, E, eqx_Du);

 // --- Submitting the terms for the governing equations ---
variable_list.set_vector_gradient_term_LHS(1,eqx_Du);

}
