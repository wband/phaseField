// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attribute_loader.h>
#include <prismspf/core/variable_container.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <typename number, typename dtType>
void
constrain_rxn_x1(const number &val, const number n, const dtType dt, number &rxn, const double &epsilon_denom)
{
  using std::max;
  using std::min;
  double epsilon = 1e-8;
  number top_lim = 1.0 - epsilon;
  number bot_lim = 0.0 + epsilon;
  number test = val + rxn * dt/(n + epsilon_denom);
  number remainder_top = max(test - top_lim, number(0.0));
  number remainder_bot = min(test + bot_lim, number(0.0));
  rxn  = (rxn * dt/(n + epsilon_denom)- remainder_top - remainder_bot) * (n + epsilon_denom) / dt;
}

template <typename number, typename dtType>
void
constrain_rxn_x2(const number &val, const number n, const dtType dt, number &rxn, const double &epsilon_denom)
{
  using std::max;
  using std::min;
  double epsilon = 1e-8;
  number top_lim = 1.0 - epsilon;
  number bot_lim = 0.0 + epsilon;
  number test = val - rxn * dt/(1.0 - n + epsilon_denom);
  number remainder_top = max(test - top_lim, number(0.0));
  number remainder_bot = min(test + bot_lim, number(0.0));
  rxn  = -(-rxn * dt/(1.0 - n + epsilon_denom) - remainder_top - remainder_bot) * (1.0 - n + epsilon_denom) / dt;
}

template <typename number, typename dtType>
void
constrain_rxn_phi(const number &val, const dtType dt, number &rxn)
{
  using std::max;
  using std::min;
  number top_lim = 1.0;
  number bot_lim = 0.0;
  number remainder_top = val + rxn * dt - top_lim;
  number remainder_bot = val + rxn * dt + bot_lim;
  remainder_top = max(remainder_top, number(0.0));
  remainder_bot = min(remainder_bot, number(0.0));
  rxn  = (rxn * dt - remainder_top - remainder_bot) / dt;
}

void
CustomAttributeLoader::load_variable_attributes()
{
  std::string dependency_string = "n, x1, x2, rxn, rxn_mu, diff1, diff2, grad(n), grad(x1), grad(x2)";
  set_variable_name(0, "n");
  set_variable_type(0, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(0, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(0, dependency_string);
  set_solve_block(0, 2);

  set_variable_name(1, "x1");
  set_variable_type(1, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(1, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(1, dependency_string);
  set_solve_block(1, 2);

  set_variable_name(2, "x2");
  set_variable_type(2, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(2, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(2, dependency_string);
  set_solve_block(2, 2);

  set_variable_name(3, "rxn");
  set_variable_type(3, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(3, Auxiliary);
  set_dependencies_value_term_rhs(3, dependency_string);
  set_solve_block(3, 1);

  set_variable_name(4, "diff1");
  set_variable_type(4, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(4, Auxiliary);
  set_dependencies_value_term_rhs(4, dependency_string);
  set_dependencies_gradient_term_rhs(4, dependency_string);
  set_solve_block(4, 0);

  set_variable_name(5, "diff2");
  set_variable_type(5, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(5, Auxiliary);
  set_dependencies_value_term_rhs(5, dependency_string);
  set_dependencies_gradient_term_rhs(5, dependency_string);
  set_solve_block(5, 0);

  set_variable_name(6, "rxn_mu");
  set_variable_type(6, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(6, Auxiliary);
  set_dependencies_value_term_rhs(6, dependency_string);
  set_dependencies_gradient_term_rhs(6, dependency_string);
  set_solve_block(6, 0);
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block) const
{
  if(solve_block == 2)
    {
      number dt = get_timestep();
      ScalarValue n  = variable_list.template get_value<ScalarValue>(0);
      ScalarValue x1  = variable_list.template get_value<ScalarValue>(1);
      ScalarValue x2  = variable_list.template get_value<ScalarValue>(2);
      ScalarValue rxn  = variable_list.template get_value<ScalarValue>(3);
      ScalarValue diff1  = variable_list.template get_value<ScalarValue>(4);
      ScalarValue diff2  = variable_list.template get_value<ScalarValue>(5);
      ScalarValue rxn_mu  = variable_list.template get_value<ScalarValue>(6);

      variable_list.set_value_term(0, n + dt * rxn);
      variable_list.set_value_term(1, x1 + dt * (diff1 + rxn));
      variable_list.set_value_term(2, x2 + dt * (diff2 - rxn));
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block,
  [[maybe_unused]] Types::Index                           index) const
{
  if(solve_block == 0)
    {
      ScalarValue n  = variable_list.template get_value<ScalarValue>(0);
      ScalarGrad  n_grad = variable_list.template get_gradient<ScalarGrad>(0);
      ScalarValue x1  = variable_list.template get_value<ScalarValue>(1);
      ScalarGrad  x1_grad = variable_list.template get_gradient<ScalarGrad>(1);
      ScalarValue x2  = variable_list.template get_value<ScalarValue>(2);
      ScalarGrad  x2_grad = variable_list.template get_gradient<ScalarGrad>(2);

      // rxn_mu
      ScalarValue rxn_mu_val = std::log(x1) - std::log(x2) + deltaG + dw_coeff * (1.0 - 2.0 * n);
      variable_list.set_value_term(6, rxn_mu_val);
      variable_list.set_gradient_term(6, n_grad * grad_coeff);
      // diff1
      variable_list.set_value_term(4, D1 * x1_grad * n_grad/(n + epsilon_denom));
      variable_list.set_gradient_term(4, -D1 * x1_grad);
      // diff2
      variable_list.set_value_term(5, -D2 * x2_grad * n_grad/(1.0 - n + epsilon_denom));
      variable_list.set_gradient_term(5, -D2 * x2_grad);
    }
  else if(solve_block == 1) // rxn
    {
      number dt = get_timestep();
      ScalarValue n  = variable_list.template get_value<ScalarValue>(0);
      ScalarValue x1  = variable_list.template get_value<ScalarValue>(1);
      ScalarValue x2  = variable_list.template get_value<ScalarValue>(2);
      ScalarValue diff1  = variable_list.template get_value<ScalarValue>(4);
      ScalarValue diff2  = variable_list.template get_value<ScalarValue>(5);
      ScalarValue rxn_mu  = variable_list.template get_value<ScalarValue>(6);

      ScalarValue rxn_val = -n * (1.0 - n) * rxn_mu;
      ScalarValue x1_temp = x1 + dt*diff1;
      ScalarValue x2_temp = x2 + dt*diff2;
      constrain_rxn_x1(x1_temp, n, dt, rxn_val, epsilon_denom);
      constrain_rxn_x2(x2_temp, n, dt, rxn_val, epsilon_denom);
      constrain_rxn_phi(n, dt, rxn_val);
      variable_list.set_value_term(3, rxn_val);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_lhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block,
  [[maybe_unused]] Types::Index                           index) const
{}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_postprocess_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block) const
{}

#include "custom_pde.inst"

PRISMS_PF_END_NAMESPACE
