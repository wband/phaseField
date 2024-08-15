// =============================================================================================
// loadPostProcessorVariableAttributes: Set the attributes of the postprocessing variables
// =============================================================================================
// This function is analogous to 'loadVariableAttributes' in 'equations.h', but for
// the postprocessing expressions. It sets the attributes for each postprocessing
// expression, including its name, whether it is a vector or scalar (only scalars are
// supported at present), its dependencies on other variables and their derivatives,
// and whether to calculate an integral of the postprocessed quantity over the entire
// domain. Note: this function is not a member of customPDE.

void variableAttributeLoader::loadPostProcessorVariableAttributes(){
	// Variable 0
	set_variable_name				(0,"f_tot");
	set_variable_type				(0,SCALAR);
        set_dependencies_value_term_RHS(0, "phi, grad(phi), grad(u)");
        set_dependencies_gradient_term_RHS(0, "");
        set_output_integral         	(0,true);

	// Variable 1
	set_variable_name				(1,"s11");
	set_variable_type				(1,SCALAR);
        set_dependencies_value_term_RHS(1, "phi, grad(u)");
        set_dependencies_gradient_term_RHS(1, "");
        set_output_integral         	(1,true);
	// Variable 2
	set_variable_name				(2,"s12");
	set_variable_type				(2,SCALAR);
        set_dependencies_value_term_RHS(2, "phi, grad(u)");
        set_dependencies_gradient_term_RHS(2, "");
        set_output_integral         	(2,true);
	// Variable 3
	set_variable_name				(3,"s22");
	set_variable_type				(3,SCALAR);
        set_dependencies_value_term_RHS(3, "phi, grad(u)");
        set_dependencies_gradient_term_RHS(3, "");
        set_output_integral         	(3,true);
	// Variable 4
	set_variable_name				(4,"e22");
	set_variable_type				(4,SCALAR);
        set_dependencies_value_term_RHS(4, "phi, grad(u)");
        set_dependencies_gradient_term_RHS(4, "");
        set_output_integral         	(4,true);
}

// =================================================================================
// Define the expressions for the post-processed fields
// =================================================================================

template <int dim,int degree>
void customPDE<dim,degree>::postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
												const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

        // --- Getting the values and derivatives of the model variables ---

		// The order parameter and its derivatives
		scalarvalueType phi = variable_list.get_scalar_value(1);
		scalargradType nx = variable_list.get_scalar_gradient(1);

		// The derivative of the displacement vector
		vectorgradType ux = variable_list.get_vector_gradient(0);
		

        // --- Setting the expressions for the terms in the postprocessing expressions ---


		// Free energy expressions and interpolation functions
                scalarvalueType f_tot = constV(0.0);


		//compute strain and stress
		dealii::VectorizedArray<double> E[dim][dim], S[dim][dim];

		for (unsigned int i=0; i<dim; i++){
		  for (unsigned int j=0; j<dim; j++){
			  E[i][j]= constV(0.5)*(ux[i][j]+ux[j][i]);

		  }
		}

for (unsigned int i=0; i<dim; i++){                     
for (unsigned int j=0; j<dim; j++){             
for (unsigned int k=0; k<dim; k++){
for (unsigned int l=0; l<dim; l++){
   S[i][j] += phi*Cijkl0[i][j][k][l]*E[k][l];
}
}
}
}


		scalarvalueType f_el = constV(0.0);

		scalarvalueType s11 = S[0][0];
		scalarvalueType s12 = S[0][1];
		scalarvalueType s22 = S[1][1];
		scalarvalueType e22 = ux[1][1];

// --- Submitting the terms for the postprocessing expressions ---

pp_variable_list.set_scalar_value_term_RHS(0, f_tot);
pp_variable_list.set_scalar_value_term_RHS(1, s11);
pp_variable_list.set_scalar_value_term_RHS(2, s12);
pp_variable_list.set_scalar_value_term_RHS(3, s22);
pp_variable_list.set_scalar_value_term_RHS(4, e22);

}
