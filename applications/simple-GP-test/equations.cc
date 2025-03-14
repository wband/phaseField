// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
// This function sets attributes for each variable/equation in the app. The
// attributes are set via standardized function calls. The first parameter for
// each function call is the variable index (starting at zero). The first set of
// variable/equation attributes are the variable name (any string), the variable
// type (SCALAR/VECTOR), and the equation type (EXPLICIT_TIME_DEPENDENT/
// TIME_INDEPENDENT/AUXILIARY). The next set of attributes describe the
// dependencies for the governing equation on the values and derivatives of the
// other variables for the value term and gradient term of the RHS and the LHS.
// The final pair of attributes determine whether a variable represents a field
// that can nucleate and whether the value of the field is needed for nucleation
// rate calculations.

void
customAttributeLoader::loadVariableAttributes()
{
  // Variable 0
  set_variable_name(0, "mu");
  set_variable_type(0, SCALAR);
  set_variable_equation_type(0, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(0, "mu, ns, nl, dnsdt, dnldt");
  set_dependencies_gradient_term_RHS(0, "grad(mu)");

  // Variable 1
  set_variable_name(1, "ns");
  set_variable_type(1, SCALAR);
  set_variable_equation_type(1, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(1, "ns, dnsdt");
  set_dependencies_gradient_term_RHS(1, "");

  // Variable 1
  set_variable_name(2, "nl");
  set_variable_type(2, SCALAR);
  set_variable_equation_type(2, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(2, "nl, dnldt");
  set_dependencies_gradient_term_RHS(2, "");

  // Variable 1
  set_variable_name(3, "dnsdt");
  set_variable_type(3, SCALAR);
  set_variable_equation_type(3, AUXILIARY);

  set_dependencies_value_term_RHS(3, "mu, ns, nl");
  set_dependencies_gradient_term_RHS(3, "grad(ns)");

  // Variable 1
  set_variable_name(4, "dnldt");
  set_variable_type(4, SCALAR);
  set_variable_equation_type(4, AUXILIARY);

  set_dependencies_value_term_RHS(4, "mu, ns, nl");
  set_dependencies_gradient_term_RHS(4, "grad(nl)");

}

// =============================================================================================
// explicitEquationRHS (needed only if one or more equation is explict time
// dependent)
// =============================================================================================
// This function calculates the right-hand-side of the explicit time-dependent
// equations for each variable. It takes "variable_list" as an input, which is a
// list of the value and derivatives of each of the variables at a specific
// quadrature point. The (x,y,z) location of that quadrature point is given by
// "q_point_loc". The function outputs two terms to variable_list -- one
// proportional to the test function and one proportional to the gradient of the
// test function. The index for each variable in this list corresponds to the
// index given at the top of this file.

template <int dim, int degree>
void
customPDE<dim, degree>::explicitEquationRHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{
  // --- Getting the values and derivatives of the model variables ---

  // c
  scalarvalueType mu  = variable_list.get_scalar_value(0);
  scalargradType  mux = variable_list.get_scalar_gradient(0);

  // n
  scalarvalueType ns  = variable_list.get_scalar_value(1);
  scalarvalueType nl  = variable_list.get_scalar_value(2);
  scalarvalueType dnsdt  = variable_list.get_scalar_value(3);
  scalarvalueType dnldt  = variable_list.get_scalar_value(4);

  // --- Setting the expressions for the terms in the governing equations ---
  scalarvalueType denom = ns*ns + nl*nl;
  scalarvalueType hs = ns*ns/denom;
  scalarvalueType hl = nl*nl/denom;

  scalarvalueType dhsdns = (2.0*ns - 2.0*ns*hs) / denom;
  scalarvalueType dhldns = ( - 2.0*ns*hl) / denom;
  scalarvalueType dhsdnl = ( - 2.0*nl*hs) / denom;
  scalarvalueType dhldnl = (2.0*nl - 2.0*nl*hl) / denom;
  
  scalarvalueType chiAA = hs*(1.0/ks) + hl*(1.0/kl);

  scalarvalueType dmudt_val = ((dhsdns*dnsdt+dhsdnl*dnldt)*(mu/ks + cs0) + 
                                (dhldns*dnsdt+dhldnl*dnldt)*(mu/kl + cl0))/chiAA;

  // Terms for the equation to evolve the concentration
  variable_list.set_scalar_value_term_RHS(0, mu + dmudt_val*userInputs.dtValue);
  variable_list.set_scalar_gradient_term_RHS(0, -mux*userInputs.dtValue);

  // Terms for the equation to evolve the order parameter
  variable_list.set_scalar_value_term_RHS(1, ns + userInputs.dtValue*dnsdt);
  variable_list.set_scalar_value_term_RHS(2, nl + userInputs.dtValue*dnldt);
}

// =============================================================================================
// nonExplicitEquationRHS (needed only if one or more equation is time
// independent or auxiliary)
// =============================================================================================
// This function calculates the right-hand-side of all of the equations that are
// not explicit time-dependent equations. It takes "variable_list" as an input,
// which is a list of the value and derivatives of each of the variables at a
// specific quadrature point. The (x,y,z) location of that quadrature point is
// given by "q_point_loc". The function outputs two terms to variable_list --
// one proportional to the test function and one proportional to the gradient of
// the test function. The index for each variable in this list corresponds to
// the index given at the top of this file.

template <int dim, int degree>
void
customPDE<dim, degree>::nonExplicitEquationRHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{
  // n
  scalarvalueType mu  = variable_list.get_scalar_value(0);
  scalarvalueType ns  = variable_list.get_scalar_value(1);
  scalarvalueType nl  = variable_list.get_scalar_value(2);
  
  scalargradType  nsx = variable_list.get_scalar_gradient(1);
  scalargradType  nlx = variable_list.get_scalar_gradient(2);
  
  scalarvalueType denom = ns*ns + nl*nl;
  scalarvalueType hs = ns*ns/denom;
  scalarvalueType hl = nl*nl/denom;

  scalarvalueType dhsdns = (2.0*ns - 2.0*ns*hs) / denom;
  scalarvalueType dhldns = ( - 2.0*ns*hl) / denom;
  scalarvalueType dhsdnl = ( - 2.0*nl*hs) / denom;
  scalarvalueType dhldnl = (2.0*nl - 2.0*nl*hl) / denom;
  
  scalarvalueType w_chems = -0.5*mu*mu/ks - mu*cs0;
  scalarvalueType w_cheml = -0.5*mu*mu/kl - mu*cl0;
  
  scalarvalueType dnsdt_val = -L*(m*(ns*ns*ns - ns + gamma*2.0*ns*nl*nl) 
      + w_chems*dhsdns + w_cheml*dhldns);
  
  scalarvalueType dnldt_val = -L*(m*(nl*nl*nl - nl + gamma*2.0*nl*ns*ns) 
      + w_chems*dhsdnl + w_cheml*dhldnl);
  
  variable_list.set_scalar_value_term_RHS(3, dnsdt_val);
  variable_list.set_scalar_value_term_RHS(4, dnldt_val);
  variable_list.set_scalar_gradient_term_RHS(3, -L*m*kappa*nsx);
  variable_list.set_scalar_gradient_term_RHS(4, -L*m*kappa*nlx);
}

// =============================================================================================
// equationLHS (needed only if at least one equation is time independent)
// =============================================================================================
// This function calculates the left-hand-side of time-independent equations. It
// takes "variable_list" as an input, which is a list of the value and
// derivatives of each of the variables at a specific quadrature point. The
// (x,y,z) location of that quadrature point is given by "q_point_loc". The
// function outputs two terms to variable_list -- one proportional to the test
// function and one proportional to the gradient of the test function -- for the
// left-hand-side of the equation. The index for each variable in this list
// corresponds to the index given at the top of this file. If there are multiple
// elliptic equations, conditional statements should be sed to ensure that the
// correct residual is being submitted. The index of the field being solved can
// be accessed by "this->currentFieldIndex".

template <int dim, int degree>
void
customPDE<dim, degree>::equationLHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{}
