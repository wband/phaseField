// =============================================================================================
// loadPostProcessorVariableAttributes: Set the attributes of the postprocessing variables
// =============================================================================================
// This function is analogous to 'loadVariableAttributes' in 'equations.h', but for
// the postprocessing expressions. It sets the attributes for each postprocessing
// expression, including its name, whether it is a vector or scalar (only scalars are
// supported at present), its dependencies on other variables and their derivatives,
// and whether to calculate an integral of the postprocessed quantity over the entire
// domain. Note: this function is not a member of customPDE.

void
customAttributeLoader::loadPostProcessorVariableAttributes()
{
	// Variable 0
	set_variable_name				(0,"f_tot");
	set_variable_type				(0,SCALAR);
        set_dependencies_value_term_RHS(0, "n, grad(n), grad(u), Ex, Gx");
        set_dependencies_gradient_term_RHS(0, "");
        set_output_integral         	(0,true);

	// Variable 1
	set_variable_name				(1,"s11");
	set_variable_type				(1,SCALAR);
        set_dependencies_value_term_RHS(1, "n, grad(u), Ex");
        set_dependencies_gradient_term_RHS(1, "");
        set_output_integral         	(1,true);
	// Variable 2
	set_variable_name				(2,"s12");
	set_variable_type				(2,SCALAR);
        set_dependencies_value_term_RHS(2, "n, grad(u), Ex");
        set_dependencies_gradient_term_RHS(2, "");
        set_output_integral         	(2,true);
	// Variable 3
	set_variable_name				(3,"s22");
	set_variable_type				(3,SCALAR);
        set_dependencies_value_term_RHS(3, "n, grad(u), Ex");
        set_dependencies_gradient_term_RHS(3, "");
        set_output_integral         	(3,true);
	// Variable 4
	set_variable_name				(4,"e22");
	set_variable_type				(4,SCALAR);
        set_dependencies_value_term_RHS(4, "grad(u), Ex");
        set_dependencies_gradient_term_RHS(4, "");
        set_output_integral         	(4,true);
	// Variable 5
	set_variable_name				(5,"f_int");
	set_variable_type				(5,SCALAR);
        set_dependencies_value_term_RHS(5, "n, grad(n), grad(u), Ex, Gx");
        set_dependencies_gradient_term_RHS(5, "");
        set_output_integral         	(5,true);
	// Variable 6
	set_variable_name				(6,"f_el");
	set_variable_type				(6,SCALAR);
        set_dependencies_value_term_RHS(6, "n, grad(u), Ex");
        set_dependencies_gradient_term_RHS(6, "");
        set_output_integral         	(6,true);

}

// =================================================================================
// Define the expressions for the post-processed fields
// =================================================================================
template <int dim, int degree>
void
customPDE<dim, degree>::postProcessedFields(
  [[maybe_unused]] const variableContainer<dim, degree, VectorizedArray<double>>
    &variable_list,
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                            &pp_variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
  [[maybe_unused]] const VectorizedArray<double>             element_volume) const
{

        // --- Getting the values and derivatives of the model variables ---

		// The order parameter and its derivatives
		scalarvalueType n = variable_list.get_scalar_value(0);
		scalargradType nx = variable_list.get_scalar_gradient(0);

		// The derivative of the displacement vector
		vectorgradType ux = variable_list.get_vector_gradient(1);
		
		// Spatially varying Young's modulus and fracture toughness
		scalarvalueType Ex = variable_list.get_scalar_value(3);
                scalarvalueType Gx = variable_list.get_scalar_value(4);
        // --- Setting the expressions for the terms in the postprocessing expressions ---
                scalarvalueType f_int = Gc0*n*Gx*3.0/8.0/ell;

		for (int i=0; i<dim; i++){
		    f_int += Gc0*Gx*0.5*ell*nx[i]*nx[i];
		}

		//compute strain and stress
		dealii::VectorizedArray<double> E[dim][dim], S[dim][dim];

		for (unsigned int i=0; i<dim; i++){
		  for (unsigned int j=0; j<dim; j++){
			  E[i][j]= 0.5*(ux[i][j]+ux[j][i]);

		  }
		}
		dealii::VectorizedArray<double> CIJ_combined[CIJ_tensor_size][CIJ_tensor_size];
		for (unsigned int i=0; i<2*dim-1+dim/3; i++){
			  for (unsigned int j=0; j<2*dim-1+dim/3; j++){
				CIJ_combined[i][j] = CIJ_base[i][j]*Ex*(1.0-2.0*n+n*n);
			}
		  }
		computeStress<dim>(CIJ_combined, E, S);
                scalarvalueType f_el = 0.0;
		for (unsigned int i=0; i<dim; i++){
		  for (unsigned int j=0; j<dim; j++){
			  f_el += 0.5 * S[i][j]*E[i][j];
		  }
		}

		scalarvalueType f_tot = f_el+f_int;
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
pp_variable_list.set_scalar_value_term_RHS(5, f_int);
pp_variable_list.set_scalar_value_term_RHS(6, f_el);

}
