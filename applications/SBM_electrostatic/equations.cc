// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attribute_loader.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

void
CustomAttributeLoader::load_variable_attributes()
{
  set_variable_name(0, "phie");
  set_variable_type(0, Scalar);
  set_variable_equation_type(0, TimeIndependent);
  set_dependencies_value_term_rhs(0, "psi,grad(psi), phie");
  set_dependencies_gradient_term_rhs(0, "psi, grad(phie)");
  set_dependencies_value_term_lhs(0, "psi,grad(psi), change(phie)");
  set_dependencies_gradient_term_lhs(0, "psi, grad(change(phie))");
  
  set_variable_name(1, "psi");
  set_variable_type(1, Scalar);
  set_variable_equation_type(1, Constant);

}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] Types::Index current_index) const
{
  
  ScalarValue phie  = variable_list.get_scalar_value(0);
	ScalarGrad  phiex = variable_list.get_scalar_gradient(0);
	ScalarValue psi  = variable_list.get_scalar_value(1);
	ScalarGrad  psix = variable_list.get_scalar_gradient(1);
	
	ScalarValue eq_phie  = psix.norm() * McV * phie * psi;
	ScalarGrad  eqx_phie =  (1.0 - psi) * KcV * phiex;
	
	variable_list.set_scalar_value_term(0, eq_phie);
	variable_list.set_scalar_gradient_term(0, eqx_phie);
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_lhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] Types::Index current_index) const
{
  ScalarValue Dphie  = variable_list.get_scalar_value(0, Change);
  ScalarGrad  Dphiex = variable_list.get_scalar_gradient(0, Change);
  ScalarValue psi  = variable_list.get_scalar_value(1);
  ScalarGrad  psix = variable_list.get_scalar_gradient(1);
  
	ScalarValue eq_phie  = -psix.norm() * McV * Dphie * psi;
	ScalarGrad eqx_phie = -(1.0 - psi) * KcV * Dphiex;

  variable_list.set_scalar_value_term(0, eq_phie, Change);
	variable_list.set_scalar_gradient_term(0,eqx_phie,Change);
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_postprocess_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{}

INSTANTIATE_TRI_TEMPLATE(CustomPDE)

PRISMS_PF_END_NAMESPACE
