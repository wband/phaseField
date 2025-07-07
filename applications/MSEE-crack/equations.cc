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
  
  set_variable_name(0, "Ex");
  set_variable_type(0, Scalar);
  set_variable_equation_type(0, Constant);

  set_variable_name(1, "u");
  set_variable_type(1, Vector);
  set_variable_equation_type(1, TimeIndependent);

  set_dependencies_value_term_rhs(1, "");
  set_dependencies_gradient_term_rhs(1, "grad(u), Ex, eigx");
  set_dependencies_value_term_lhs(1, "");
  set_dependencies_gradient_term_lhs(1, "grad(change(u)), Ex");
  
  set_variable_name(2, "eigx");
  set_variable_type(2, Scalar);
  set_variable_equation_type(2, Constant);

  set_variable_name(3, "s11");
  set_variable_type(3, Scalar);
  set_variable_equation_type(3, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(3, "Ex, grad(u), eigx");
  set_dependencies_gradient_term_rhs(3, "");
  set_is_postprocessed_field(3, true);

  set_variable_name(4, "s12");
  set_variable_type(4, Scalar);
  set_variable_equation_type(4, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(4, "Ex, grad(u), eigx");
  set_dependencies_gradient_term_rhs(4, "");
  set_is_postprocessed_field(4, true);

  set_variable_name(5, "s22");
  set_variable_type(5, Scalar);
  set_variable_equation_type(5, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(5, "Ex, grad(u), eigx");
  set_dependencies_gradient_term_rhs(5, "");
  set_is_postprocessed_field(5, true);

}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] Types::Index current_index) const
{
  if (current_index == 1)
    {
      ScalarValue Ex = variable_list.get_scalar_value(0);
      ScalarValue eigx = variable_list.get_scalar_value(2);
      VectorGrad ux = variable_list.get_vector_symmetric_gradient(1);
      for (unsigned int i = 0; i < dim; i++)
        {
          ux[i][i] -= eigx;
        }
      VectorGrad stress;
      compute_stress<dim, ScalarValue>(stiffness, Ex*ux, stress);
      variable_list.set_vector_gradient_term(1, stress);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_lhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] Types::Index current_index) const
{
  if (current_index == 1)
    {
      ScalarValue Ex = variable_list.get_scalar_value(0);
      VectorGrad change_ux = variable_list.get_vector_symmetric_gradient(1, Change);
      VectorGrad stress;
      compute_stress<dim, ScalarValue>(stiffness, Ex*change_ux, stress);
      variable_list.set_vector_gradient_term(1, -stress, Change);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_postprocess_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{
  ScalarValue Ex  = variable_list.get_scalar_value(0);
  VectorGrad  ux = variable_list.get_vector_symmetric_gradient(1);
  ScalarValue eigx = variable_list.get_scalar_value(2);
  VectorGrad stress;
  
  for (unsigned int i = 0; i < dim; i++)
    {
      ux[i][i] -= eigx;
    }
  compute_stress<dim, ScalarValue>(stiffness, Ex*ux, stress);
  
  ScalarValue s11   = stress[0][0];
  ScalarValue s12   = stress[0][1];
  ScalarValue s22   = stress[1][1];

  variable_list.set_scalar_value_term(3, s11);
  variable_list.set_scalar_value_term(4, s12);
  variable_list.set_scalar_value_term(5, s22);
}

INSTANTIATE_TRI_TEMPLATE(CustomPDE)

PRISMS_PF_END_NAMESPACE
