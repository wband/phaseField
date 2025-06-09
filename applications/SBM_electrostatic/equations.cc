// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attribute_loader.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

void
customAttributeLoader::loadVariableAttributes()
{
  set_variable_name(0, "phie");
  set_variable_type(0, SCALAR);
  set_variable_equation_type(0, TIME_INDEPENDENT);
  set_dependencies_value_term_RHS(0, "psi,grad(psi), phie");
  set_dependencies_gradient_term_RHS(0, "psi, grad(phie)");
  set_dependencies_value_term_LHS(0, "psi,grad(psi), change(phie)");
  set_dependencies_gradient_term_LHS(0, "psi, grad(change(phie))");
  
  set_variable_name(1, "psi");
  set_variable_type(1, SCALAR);
  set_variable_equation_type(1, CONSTANT);

}

template <unsigned int dim, unsigned int degree, typename number>
void
customPDE<dim, degree, number>::compute_explicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{
}

template <unsigned int dim, unsigned int degree, typename number>
void
customPDE<dim, degree, number>::compute_nonexplicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] types::index current_index) const
{
  
  scalarValue phie  = variable_list.get_scalar_value(0);
	scalarGrad  phiex = variable_list.get_scalar_gradient(0);
	scalarValue psi  = variable_list.get_scalar_value(1);
	scalarGrad  psix = variable_list.get_scalar_gradient(1);
	
	scalarValue eq_phie  = psix.norm() * McV * phie * psi;
	scalarGrad  eqx_phie =  (1.0 - psi) * KcV * phiex;
	
	variable_list.set_scalar_value_term(0, eq_phie);
	variable_list.set_scalar_gradient_term(0, eqx_phie);
}

template <unsigned int dim, unsigned int degree, typename number>
void
customPDE<dim, degree, number>::compute_nonexplicit_LHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] types::index current_index) const
{
  scalarValue Dphie  = variable_list.get_scalar_value(0, CHANGE);
  scalarGrad  Dphiex = variable_list.get_scalar_gradient(0, CHANGE);
  scalarValue psi  = variable_list.get_scalar_value(1);
  scalarGrad  psix = variable_list.get_scalar_gradient(1);
  
	scalarValue eq_phie  = -psix.norm() * McV * Dphie * psi;
	scalarGrad eqx_phie = -(1.0 - psi) * KcV * Dphiex;

  variable_list.set_scalar_value_term(0, eq_phie, CHANGE);
	variable_list.set_scalar_gradient_term(0,eqx_phie,CHANGE);
}

template <unsigned int dim, unsigned int degree, typename number>
void
customPDE<dim, degree, number>::compute_postprocess_explicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{}

INSTANTIATE_TRI_TEMPLATE(customPDE)

PRISMS_PF_END_NAMESPACE
