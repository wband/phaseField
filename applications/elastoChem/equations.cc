// List of variables and residual equations for the mechanics example application

// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
void variableAttributeLoader::loadVariableAttributes(){

// Variable 0 - concentration field
	set_variable_name				(0,"n");
	set_variable_type				(0,SCALAR);
	set_variable_equation_type		(0,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(0, "n, grad(u), phi");
    set_dependencies_gradient_term_RHS(0, "");
	
// Variable 1 - displacement field
	set_variable_name				(1,"u");
	set_variable_type				(1,VECTOR);
	set_variable_equation_type		(1,TIME_INDEPENDENT);

    set_dependencies_value_term_RHS(1, "");
    set_dependencies_gradient_term_RHS(1, "n, grad(u), phi");
    set_dependencies_value_term_LHS(1, "");
    set_dependencies_gradient_term_LHS(1, "n, grad(change(u)), phi");

// Variable 2 - phase field
	set_variable_name				(2,"phi");
	set_variable_type				(2,SCALAR);
	set_variable_equation_type		(2,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(2, "phi");
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
// The phase field and the variable containing its update
scalarvalueType phi = variable_list.get_scalar_value(2);
//n
scalarvalueType n = variable_list.get_scalar_value(0);
//u
vectorgradType ux = variable_list.get_vector_gradient(1);

//compute strain tensor
vectorgradType nonchem_strain, stress;
 for (unsigned int i=0; i<dim; i++){
 for (unsigned int j=0; j<dim; j++){
          nonchem_strain[i][j]= 0.5*(ux[i][j]+ux[j][i]);
    if (i == j){
        nonchem_strain[i][j] -= coupling*n;
    }
 }
 }

scalarvalueType elastic_potential;

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
for (unsigned int k=0; k<dim; k++){
for (unsigned int l=0; l<dim; l++){
   stress[i][j] += (phi*Cijkl1[i][j][k][l] + (1.0-phi)*Cijkl0[i][j][k][l])*
       nonchem_strain[k][l];
   if (k == l){
       elastic_potential += stress[i][j]*coupling;
   }
}
}
}
}


variable_list.set_scalar_value_term_RHS(0, n-userInputs.dtValue*(kwell*n - elastic_potential - mu_bar));
variable_list.set_scalar_value_term_RHS(2,phi);
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
scalarvalueType phi = variable_list.get_scalar_value(2);
//u
vectorgradType ux = variable_list.get_vector_gradient(1);


//compute strain tensor
vectorgradType nonchem_strain, stress;
 for (unsigned int i=0; i<dim; i++){
 for (unsigned int j=0; j<dim; j++){
          nonchem_strain[i][j]= 0.5*(ux[i][j]+ux[j][i]);
    if (i == j){
        nonchem_strain[i][j] -= coupling*n;
    }
 }
 }

scalarvalueType elastic_potential; 

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
for (unsigned int k=0; k<dim; k++){
for (unsigned int l=0; l<dim; l++){
   stress[i][j] += (phi*Cijkl1[i][j][k][l] + (1.0-phi)*Cijkl0[i][j][k][l])*
       nonchem_strain[k][l];
   if (k == l){
       elastic_potential += stress[i][j]*coupling;
   }
}
}
}
}

// implied: f = 0.5*k*n**2 + 0.5*C*(e-ec*n)**2
// mu = k*n - C*e*ec + C*ec*ec*n = k*n - C*(e-ec*n)*ec
// set mu constant: mu_bar
// n = (mu_bar + elastic_potential)/k

// --- Submitting the terms for the governing equations ---
// phase field
//variable_list.set_scalar_value_term_RHS(0, n-userInputs.dtValue*(kwell*n - elastic_potential - mu_bar));
// mechanics
variable_list.set_vector_gradient_term_RHS(1,stress);

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

// --- Setting the expressions for the terms in the governing equations ---
//n
//scalarvalueType Dn = variable_list.get_change_in_scalar_value(0);
scalarvalueType Dn = variable_list.get_scalar_value(0);
//u
vectorgradType Dux = variable_list.get_change_in_vector_gradient(1);

scalarvalueType phi = variable_list.get_scalar_value(2);

// mu_bar = k*n - C*(e-ec*n)*ec
// 0 = -(k*dn + C*ec*ec*dn - C*ec*de)
// 0 = -(k*dn + C*ec*(ec*dn - de))
// n = (mu_bar + elastic_potential)/k
// C*(e-ec*n)
// C*(de-ec*dn)


//compute strain tensor
vectorgradType nonchem_strain, stress;
 for (unsigned int i=0; i<dim; i++){
 for (unsigned int j=0; j<dim; j++){
          nonchem_strain[i][j]= 0.5*(Dux[i][j]+Dux[j][i]);
    if (i == j){
        nonchem_strain[i][j] -= coupling*Dn;
    }
 }
 }

scalarvalueType elastic_potential; 

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
for (unsigned int k=0; k<dim; k++){
for (unsigned int l=0; l<dim; l++){
   stress[i][j] += (phi*Cijkl1[i][j][k][l] + (1.0-phi)*Cijkl0[i][j][k][l])*
        nonchem_strain[k][l];
//   if (k == l){
//       elastic_potential += stress[i][j]*coupling;
//}
}
}
}
}

 // --- Submitting the terms for the governing equations ---
//variable_list.set_scalar_value_term_LHS(0,-(kwell*Dn - elastic_potential));
//variable_list.set_scalar_value_term_LHS(0,-(kwell*Dn));
variable_list.set_vector_gradient_term_LHS(1,-stress);

}
