// List of variables and residual equations for the mechanics example application

// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadVariableAttributes(){

	// Variable 2
	set_variable_name				(0,"c");
	set_variable_type				(0,SCALAR);
	set_variable_equation_type		(0,TIME_INDEPENDENT);
	
	set_variable_name				(1,"phi1");
	set_variable_type				(1,SCALAR);
	set_variable_equation_type		(1, EXPLICIT_TIME_DEPENDENT);

	set_variable_name				(2,"phi2");
	set_variable_type				(2,SCALAR);
	set_variable_equation_type		(2, EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(0, "");
    set_dependencies_gradient_term_RHS(0, "grad(c), phi1, grad(phi1), phi2, grad(phi2)");
    set_dependencies_value_term_LHS(0, "");
    set_dependencies_gradient_term_LHS(0, "grad(change(c)), phi1, grad(phi1), phi2, grad(phi2)");

    set_dependencies_value_term_RHS(1, "phi1");
    set_dependencies_gradient_term_RHS(1, "");
    set_dependencies_value_term_RHS(2, "phi2");
    set_dependencies_gradient_term_RHS(2, "");
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
scalarvalueType phi1 = variable_list.get_scalar_value(1);
scalarvalueType phi2 = variable_list.get_scalar_value(2);
variable_list.set_scalar_value_term_RHS(1,phi1);
variable_list.set_scalar_value_term_RHS(2,phi2);
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
scalargradType cx = variable_list.get_scalar_gradient(0);
scalarvalueType phi1 = variable_list.get_scalar_value(1);
scalarvalueType phi2 = variable_list.get_scalar_value(2);
scalargradType phi1x = variable_list.get_scalar_gradient(1);
scalargradType phi2x = variable_list.get_scalar_gradient(2);

scalarvalueType normgrad1 = std::sqrt(phi1x.norm_square());
scalarvalueType normgrad2 = std::sqrt(phi2x.norm_square());
scalargradType n1 = phi1x/(normgrad1+constV(1.0e-16));
scalargradType n2 = phi2x/(normgrad2+constV(1.0e-16));
scalargradType Dcx;
for (unsigned int i=0; i<dim; i++){
    for (unsigned int j=0; j<dim; j++){
        dealii::VectorizedArray<double> D1,D2;
        if (i == j){
            D1 = 1.0;
            D2 = 1.0;
            D1 += -parallel*(1.0-n1[i]*n1[j])*4.0*phi1*phi2;
            D2 += -parallel*(1.0-n2[i]*n2[j])*4.0*phi1*phi2;
        } else {
            D1 = 0.0;
            D2 = 0.0;
            D1 += -parallel*(-n1[i]*n1[j])*4.0*phi1*phi2;
            D2 += -parallel*(-n2[i]*n2[j])*4.0*phi1*phi2;
        }
        D1 += -orthogonal*n1[i]*n1[j]*4.0*phi1*phi2;
        D2 += -orthogonal*n2[i]*n2[j]*4.0*phi1*phi2;
        Dcx[i] += (D1*phi1 + D2*phi2)*cx[j];
    }
}
variable_list.set_scalar_gradient_term_RHS(0,-Dcx);
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
//u
scalargradType cx = variable_list.get_change_in_scalar_gradient(0);
scalarvalueType phi1 = variable_list.get_scalar_value(1);
scalarvalueType phi2 = variable_list.get_scalar_value(2);
scalargradType phi1x = variable_list.get_scalar_gradient(1);
scalargradType phi2x = variable_list.get_scalar_gradient(2);
scalarvalueType normgrad1 = std::sqrt(phi1x.norm_square());
scalarvalueType normgrad2 = std::sqrt(phi2x.norm_square());
scalargradType n1 = phi1x/(normgrad1+constV(1.0e-16));
scalargradType n2 = phi2x/(normgrad2+constV(1.0e-16));
scalargradType Dcx;
for (unsigned int i=0; i<dim; i++){
    for (unsigned int j=0; j<dim; j++){
        dealii::VectorizedArray<double> D1,D2;
        if (i == j){
            D1 = 1.0;
            D2 = 1.0;
            D1 += -parallel*(1.0-n1[i]*n1[j])*4.0*phi1*phi2;
            D2 += -parallel*(1.0-n2[i]*n2[j])*4.0*phi1*phi2;
        } else {
            D1 = 0.0;
            D2 = 0.0;
            D1 += -parallel*(-n1[i]*n1[j])*4.0*phi1*phi2;
            D2 += -parallel*(-n2[i]*n2[j])*4.0*phi1*phi2;
        }
        D1 += -orthogonal*n1[i]*n1[j]*4.0*phi1*phi2;
        D2 += -orthogonal*n2[i]*n2[j]*4.0*phi1*phi2;
        Dcx[i] += (D1*phi1 + D2*phi2)*cx[j];
    }
}
// --- Submitting the terms for the governing equations ---
variable_list.set_scalar_gradient_term_LHS(0,Dcx);

}
