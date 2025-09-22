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
  set_variable_name(0, "c_new");
  set_variable_type(0, Scalar);
  set_variable_equation_type(0, TimeIndependent);
  set_dependencies_value_term_rhs(0, "c_new, c_old1, c_old2, grad(n), n");
  set_dependencies_gradient_term_rhs(0, "grad(c_new)");
  set_dependencies_value_term_lhs(0, "change(c_new), grad(change(c_new)), grad(n), n");
  set_dependencies_gradient_term_lhs(0, "grad(change(c_new))");

  set_variable_name(1, "c_old1");
  set_variable_type(1, Scalar);
  set_variable_equation_type(1, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(1, "c_new");
  set_dependencies_gradient_term_rhs(1, "");
  set_solve_block(1, 1);

  set_variable_name(2, "c_old2");
  set_variable_type(2, Scalar);
  set_variable_equation_type(2, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(2, "c_old1");
  set_dependencies_gradient_term_rhs(2, "");
  set_solve_block(2, 1);

  set_variable_name(3, "n");
  set_variable_type(3, Scalar);
  set_variable_equation_type(3, Constant);

}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block) const
{
  if (solve_block == 1)
    {
      variable_list.set_value_term(2, 
          variable_list.template get_value<ScalarValue>(1));
      variable_list.set_value_term(1,
          variable_list.template get_value<ScalarValue>(0));
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block,
  [[maybe_unused]] Types::Index                           current_index) const
{
  if ((current_index == 0) && (solve_block == 0))
    {
      if (this->get_user_inputs().get_temporal_discretization().get_time() == 0.0)
        {
          // bootstrap with one iteration of backward Euler
          ScalarValue c      = variable_list.template get_value<ScalarValue>(0);
          ScalarGrad  cx     = variable_list.template get_gradient<ScalarGrad>(0);
          ScalarValue c_old1 = variable_list.template get_value<ScalarValue>(1);
          ScalarValue n  = variable_list.template get_value<ScalarValue>(3);
          ScalarGrad  nx = variable_list.template get_gradient<ScalarGrad>(3);
          ScalarValue bc_indicator = 0.5*(1.0+std::tanh((q_point_loc[0]-channel_depth)/ell));
          ScalarValue eq_c  = c_old1 - c
              + this->get_timestep() * McV * (-c + bc_indicator*0.01)*nx.norm_square()/(n*n);
          ScalarGrad  eqx_c = -this->get_timestep() * McV * cx;
          variable_list.set_value_term(0, eq_c);
          variable_list.set_gradient_term(0, eqx_c);
        }
      else
        {
          // use BDF2 for all subsequent time steps
          ScalarValue c     = variable_list.template get_value<ScalarValue>(0);
          ScalarGrad  cx    = variable_list.template get_gradient<ScalarGrad>(0);
          ScalarValue c_old1 = variable_list.template get_value<ScalarValue>(1);
          ScalarValue c_old2 = variable_list.template get_value<ScalarValue>(2);
          ScalarValue n  = variable_list.template get_value<ScalarValue>(3);
          ScalarGrad  nx = variable_list.template get_gradient<ScalarGrad>(3);
          ScalarValue bc_indicator = 0.5*(1.0+std::tanh((q_point_loc[0]-channel_depth)/ell));
          ScalarValue eq_c  = -1.0/3.0*c_old2 + 4.0/3.0*c_old1 - c
              + 2.0/3.0*this->get_timestep() * McV * (-c + bc_indicator*0.01)*nx.norm_square()/(n*n);
          ScalarGrad  eqx_c = -2.0/3.0*this->get_timestep() * McV * cx;
          variable_list.set_value_term(0, eq_c);
          variable_list.set_gradient_term(0, eqx_c);
        }
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

  if ((current_index == 0) && (solve_block == 0))
    {
      double prefactor = 1.0;
      if (this->get_user_inputs().get_temporal_discretization().get_time() == 0.0) prefactor = 2.0/3.0;
      ScalarValue change_c  = variable_list.template get_value<ScalarValue>(0, Change);
      ScalarGrad  change_cx = variable_list.template get_gradient<ScalarGrad>(0, Change);
      ScalarValue n = variable_list.template get_value<ScalarValue>(3);
      ScalarGrad  nx    = variable_list.template get_gradient<ScalarGrad>(3);
      ScalarValue eq_change_c  = change_c*(1.0 + prefactor*this->get_timestep() * McV*nx.norm_square()/(n*n));
      ScalarGrad  eqx_change_c = prefactor*this->get_timestep() * McV * change_cx;

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
}

#include "custom_pde.inst"

PRISMS_PF_END_NAMESPACE
