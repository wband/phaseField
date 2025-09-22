// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attribute_loader.h>
#include <prismspf/core/variable_container.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

void
CustomAttributeLoader::load_variable_attributes()
{
  set_variable_name(0, "c");
  set_variable_type(0, Scalar);
  set_variable_equation_type(0, ImplicitTimeDependent);
  set_dependencies_value_term_rhs(0, "c, old_1(c), grad(n), n");
  set_dependencies_gradient_term_rhs(0, "grad(c)");
  set_dependencies_value_term_lhs(0, "change(c), grad(change(c)), grad(n), n");
  set_dependencies_gradient_term_lhs(0, "grad(change(c))");

  set_variable_name(2, "n");
  set_variable_type(2, Scalar);
  set_variable_equation_type(2, Constant);

  set_variable_name(1, "u");
  set_variable_type(1, Vector);
  set_variable_equation_type(1, TimeIndependent);
  set_dependencies_value_term_rhs(1, "");
  set_dependencies_gradient_term_rhs(1, "grad(u), n, c");
  set_dependencies_value_term_lhs(1, "");
  set_dependencies_gradient_term_lhs(1, "grad(change(u)), n");

  set_variable_name(3, "s11");
  set_variable_type(3, Scalar);
  set_variable_equation_type(3, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(3, "n, grad(u), c");
  set_dependencies_gradient_term_rhs(3, "");
  set_is_postprocessed_field(3, true);

  set_variable_name(4, "s12");
  set_variable_type(4, Scalar);
  set_variable_equation_type(4, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(4, "n, grad(u), c");
  set_dependencies_gradient_term_rhs(4, "");
  set_is_postprocessed_field(4, true);

  set_variable_name(5, "s22");
  set_variable_type(5, Scalar);
  set_variable_equation_type(5, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(5, "n, grad(u), c");
  set_dependencies_gradient_term_rhs(5, "");
  set_is_postprocessed_field(5, true);

}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block) const
{}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block,
  [[maybe_unused]] Types::Index                           current_index) const
{
  if (current_index == 1)
    {
      ScalarValue n = variable_list.template get_value<ScalarValue>(2);
      ScalarValue c = variable_list.template get_value<ScalarValue>(0);
      VectorGrad ux = variable_list.template get_symmetric_gradient<VectorGrad>(1);
      for (unsigned int i = 0; i < dim; i++)
        {
          ux[i][i] -= 0.01 * c;
        }
      VectorGrad stress;
      compute_stress<dim, ScalarValue>(stiffness, n*ux, stress);
      variable_list.set_gradient_term(1, stress);
    }
  if (current_index == 0)
    {
      ScalarValue c     = variable_list.template get_value<ScalarValue>(0);
      ScalarValue old_c = variable_list.template get_value<ScalarValue>(0, OldOne);
      ScalarGrad  cx    = variable_list.template get_gradient<ScalarGrad>(0);
      ScalarGrad  nx    = variable_list.template get_gradient<ScalarGrad>(2);
      ScalarValue n = variable_list.template get_value<ScalarValue>(2);
      ScalarValue bc_indicator = 0.5*(1.0+std::tanh((q_point_loc[0]-channel_depth)/ell));
      ScalarValue eq_c  = old_c - c
          + this->get_timestep() * McV * (-c + bc_indicator*0.01)*nx.norm_square()/(n*n);
      ScalarGrad  eqx_c = -this->get_timestep() * McV * cx;

      variable_list.set_value_term(0, eq_c);
      variable_list.set_gradient_term(0, eqx_c);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_lhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block,
  [[maybe_unused]] Types::Index                           current_index) const
{
  if (current_index == 1)
    {
      ScalarValue n = variable_list.template get_value<ScalarValue>(2);
      VectorGrad change_ux = variable_list.template get_symmetric_gradient<VectorGrad>(1,Change);
      VectorGrad stress;
      compute_stress<dim, ScalarValue>(stiffness, n*change_ux, stress);
      variable_list.set_gradient_term(1, stress, Change);
    }
  if (current_index == 0)
    {
      ScalarValue change_c  = variable_list.template get_value<ScalarValue>(0, Change);
      ScalarGrad  change_cx = variable_list.template get_gradient<ScalarGrad>(0, Change);
      ScalarGrad  nx    = variable_list.template get_gradient<ScalarGrad>(2);
      ScalarValue n = variable_list.template get_value<ScalarValue>(2);
      ScalarValue eq_change_c  = change_c*(1.0 + this->get_timestep() * McV*nx.norm_square()/(n*n));
      ScalarGrad  eqx_change_c = this->get_timestep() * McV * change_cx;

      variable_list.set_value_term(0, eq_change_c, Change);
      variable_list.set_gradient_term(0, eqx_change_c, Change);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_postprocess_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block) const
{
  ScalarValue n  = variable_list.template get_value<ScalarValue>(2);
  VectorGrad  ux = variable_list.template get_symmetric_gradient<VectorGrad>(1);
  ScalarValue c = variable_list.template get_value<ScalarValue>(0);
  VectorGrad stress;
  
  for (unsigned int i = 0; i < dim; i++)
    {
      ux[i][i] -= c;
    }
  compute_stress<dim, ScalarValue>(stiffness, n*ux, stress);
  
  ScalarValue s11   = stress[0][0];
  ScalarValue s12   = stress[0][1];
  ScalarValue s22   = stress[1][1];

  variable_list.set_value_term(3, s11);
  variable_list.set_value_term(4, s12);
  variable_list.set_value_term(5, s22);
}

#include "custom_pde.inst"

PRISMS_PF_END_NAMESPACE
