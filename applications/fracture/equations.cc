// List of variables and residual equations for the mechanics example application

// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void
customAttributeLoader::loadVariableAttributes()
{
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
    set_dependencies_gradient_term_RHS(1, "n, grad(u), Ex");
    set_dependencies_value_term_LHS(1, "");
    set_dependencies_gradient_term_LHS(1, "n, grad(change(u)), Ex");

// Variable 2 - phase field update
	set_variable_name				(2,"dndt");
	set_variable_type				(2,SCALAR);
	set_variable_equation_type		(2,AUXILIARY);

    set_dependencies_value_term_RHS(2, "n, grad(u), Ex, Gx");
    set_dependencies_gradient_term_RHS(2, "grad(n), Gx");

// Variable 3 - spatial variation in Young's modulus
	set_variable_name				(3,"Ex");
	set_variable_type				(3,SCALAR);
	set_variable_equation_type		(3,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(3, "Ex");
    set_dependencies_gradient_term_RHS(3, "");

// Variable 4 - spatial variation in fracture toughness
	set_variable_name				(4,"Gx");
	set_variable_type				(4,SCALAR);
	set_variable_equation_type		(4,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(4, "Gx");
    set_dependencies_gradient_term_RHS(4, "");
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
void
customPDE<dim, degree>::explicitEquationRHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{
// The phase field and the variable containing its update
scalarvalueType n = variable_list.get_scalar_value(0);
scalarvalueType dndt = variable_list.get_scalar_value(2);

// enforcing irreversibility: phase field must increase monotonically
for (unsigned int j=0; j<dndt.size(); ++j){
    if (dndt[j] > 0.0) dndt[j] = 0.0;
// phase field also cannot exceed 1
    if (n[j] - dndt[j]*(userInputs.dtValue) > 1.0) dndt[j] = (n[j]-1.0)/(userInputs.dtValue);
}

// update to phase field 
scalarvalueType eq_n = (n-constV(userInputs.dtValue*Mn)*dndt);

variable_list.set_scalar_value_term_RHS(0,eq_n);
variable_list.set_scalar_value_term_RHS(3, variable_list.get_scalar_value(3) );
variable_list.set_scalar_value_term_RHS(4, variable_list.get_scalar_value(4) );
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
void
customPDE<dim, degree>::nonExplicitEquationRHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{
 // --- Getting the values and derivatives of the model variables ---
//n
scalarvalueType n = variable_list.get_scalar_value(0);
scalargradType nx = variable_list.get_scalar_gradient(0);
//u
vectorgradType ux = variable_list.get_vector_gradient(1);
//material property fields
scalarvalueType Ex = variable_list.get_scalar_value(3);
scalarvalueType Gx = variable_list.get_scalar_value(4);

dealii::Tensor<2, CIJ_tensor_size, dealii::VectorizedArray<double> > CIJ;
dealii::Tensor<2, dim, dealii::VectorizedArray<double> > E;
vectorgradType S;

for (unsigned int i=0; i<2*dim-1+dim/3; i++){
          for (unsigned int j=0; j<2*dim-1+dim/3; j++){
                   CIJ[i][j] = CIJ_base[i][j]*Ex;
          }
}
E = symmetrize(ux);
//compute stress tensor
computeStress<dim>(CIJ, E, S);

//compute RHS for mechanics solver (multiply stress by interpolation function)
vectorgradType eqx_u;
for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
		eqx_u[i][j] = -S[i][j]*(1.0-2.0*n+n*n);
	}
}

// compute elastic energy density (without interpolation function)
dealii::VectorizedArray<double> elastic_energy = 0.0;

for (unsigned int i=0; i<dim; i++){
          for (unsigned int j=0; j<dim; j++){
                  elastic_energy += 0.5*S[i][j]*E[i][j];
          }
}

// phase field update: (interpolation function)*(elastic energy)+bulk energy
scalarvalueType eq_n = (2.0*(n-1.0)*elastic_energy + Gc0*Gx*3.0/8.0/ell)*Mn;

// phase field update: gradient term (Laplacian in strong form)
scalargradType eqx_n = ell*nx*Gc0*Gx*3.0/8.0*Mn;

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
void
customPDE<dim, degree>::equationLHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{
// --- Getting the values and derivatives of the model variables ---
// phase field
scalarvalueType n = variable_list.get_scalar_value(0);
// displacement
vectorgradType Dux = variable_list.get_change_in_vector_gradient(1);
// spatial variation in Young's modulus
scalarvalueType Ex = variable_list.get_scalar_value(3);

// --- Setting the expressions for the terms in the governing equations ---
vectorgradType eqx_Du;

// compute stiffness tensor
dealii::Tensor<2, CIJ_tensor_size, dealii::VectorizedArray<double> > CIJ;
for (unsigned int i=0; i<2*dim-1+dim/3; i++){
          for (unsigned int j=0; j<2*dim-1+dim/3; j++){
                  CIJ[i][j] = CIJ_base[i][j]*Ex*(1.0-2.0*n+n*n);
          }
}

//compute strain tensor
dealii::Tensor<2, dim, dealii::VectorizedArray<double> > E;
E = symmetrize(Dux);
//compute stress tensor
computeStress<dim>(CIJ, E, eqx_Du);

 // --- Submitting the terms for the governing equations ---
variable_list.set_vector_gradient_term_LHS(1,eqx_Du);

}
