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
  std::string dependency_string = "n, x1, x2, rxn, rxn_mu, grad(n), grad(x1), grad(x2)";
  set_variable_name(0, "n");
  set_variable_type(0, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(0, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(0, dependency_string);
  set_solve_block(0, 2);

  set_variable_name(1, "x1");
  set_variable_type(1, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(1, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(1, dependency_string);
  set_dependencies_gradient_term_rhs(1, dependency_string);
  set_solve_block(1, 2);

  set_variable_name(2, "x2");
  set_variable_type(2, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(2, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(2, dependency_string);
//  set_dependencies_gradient_term_rhs(2, dependency_string);
  set_solve_block(2, 2);

  set_variable_name(3, "rxn");
  set_variable_type(3, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(3, Auxiliary);
  set_dependencies_value_term_rhs(3, dependency_string);
  set_solve_block(3, 1);

  set_variable_name(4, "rxn_mu");
  set_variable_type(4, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(4, Auxiliary);
  set_dependencies_value_term_rhs(4, dependency_string);
  set_dependencies_gradient_term_rhs(4, dependency_string);
  set_solve_block(4, 0);
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
      ScalarGrad  n_grad = variable_list.template get_gradient<ScalarGrad>(0);
      ScalarValue x1  = variable_list.template get_value<ScalarValue>(1);
      ScalarGrad  x1_grad = variable_list.template get_gradient<ScalarGrad>(1);
      ScalarValue x2  = variable_list.template get_value<ScalarValue>(2);
      ScalarGrad  x2_grad = variable_list.template get_gradient<ScalarGrad>(2);
      ScalarValue rxn  = variable_list.template get_value<ScalarValue>(3);

      variable_list.set_value_term(0, n + dt * rxn);
      variable_list.set_value_term(1, x1 + dt * (D1 * x1_grad * n_grad/(n + epsilon_denom) + rxn));
      variable_list.set_gradient_term(1, dt * (-D1 * x1_grad));
      variable_list.set_value_term(2, x2);
  //    variable_list.set_gradient_term(2, ScalarGrad(0.0));
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
  using std::max;
  using std::min;
  const number upper(1.0 - 1e-4);
  const number lower(1e-4);
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
      variable_list.set_value_term(4, rxn_mu_val);
      variable_list.set_gradient_term(4, n_grad * grad_coeff);
    }
  else if(solve_block == 1) // rxn
    {
      number dt = get_timestep();
      ScalarValue n  = variable_list.template get_value<ScalarValue>(0);
      ScalarGrad  n_grad = variable_list.template get_gradient<ScalarGrad>(0);
      ScalarValue x1  = variable_list.template get_value<ScalarValue>(1);
      ScalarValue x2  = variable_list.template get_value<ScalarValue>(2);
      ScalarValue rxn_mu  = variable_list.template get_value<ScalarValue>(4);
      ScalarValue rxn_val = - n_grad.norm_square() * rxn_mu;
      ScalarValue remainder_upper;
      ScalarValue remainder_lower;
      ScalarValue test;
      // find the excess reaction that would put n, x1, and x2 out of bounds and
      // solve for the remaining reaction that would be allowed within the bounds
      // constraining for n
      test = n + rxn_val * dt;
      remainder_upper = max(test - upper, ScalarValue(0.0));
      remainder_lower = min(test - lower, ScalarValue(0.0));
      // new = n + dt*rxn => rxn = (new - n)/dt
      // for constrained cases, new = test - remainder, otherwise new = test
      rxn_val = (rxn_val * dt - remainder_upper - remainder_lower) / dt;
      // constraining for x1
 /*     test = x1 + dt * (diff1 + rxn_val);
      remainder_upper = max(test - upper, ScalarValue(0.0));
      remainder_lower = min(test - lower, ScalarValue(0.0));
      // new = x1 + dt*(diff1 + rxn) => rxn = (new - x1 - dt*diff1)/dt
      rxn_test = (rxn_val * dt - remainder_upper - remainder_lower) / dt;
      // always using the smaller magnitude reaction (looks weird bc no sign function)
      rxn_val = max(min(rxn_val,ScalarValue(0.0)),min(rxn_test,ScalarValue(0.0)))
              + min(max(rxn_val,ScalarValue(0.0)),max(rxn_test,ScalarValue(0.0))); 
      // constraining for x2
      test = x2 + dt * (diff2 - rxn_val);
      remainder_upper = max(test - upper, ScalarValue(0.0));
      remainder_lower = min(test - lower, ScalarValue(0.0));
      // new = x2 + dt*(diff2 - rxn) => rxn = -(new - x2 - dt*diff2)/dt
      rxn_test = -(-rxn_val * dt - remainder_upper - remainder_lower) / dt;
      rxn_val = max(min(rxn_val,ScalarValue(0.0)),min(rxn_test,ScalarValue(0.0)))
              + min(max(rxn_val,ScalarValue(0.0)),max(rxn_test,ScalarValue(0.0))); */

/* readable but possibly not efficient version of constraint logic */
/*      ScalarValue n_test = n + rxn_val * dt;
      for (unsigned int i = 0; i < rxn_val.size(); ++i)
        {
          if (n_test[i] > upper) rxn_val[i] = (upper-n[i])/dt;
          else if (n_test[i] < lower) rxn_val[i] = (lower-n[i])/dt;
        }
      ScalarValue x1_test = x1 + dt * (diff1 + rxn_val);
      for (unsigned int i = 0; i < rxn_val.size(); ++i)
        {
          if (x1_test[i] > upper) rxn_val[i] = (upper-x1[i]-dt*diff1[i])/dt;
          else if (x1_test[i] < lower) rxn_val[i] = (lower-x1_test[i]-dt*diff1[i])/dt;
        }
      ScalarValue x2_test = x2 + dt * (diff2 - rxn_val);
      for (unsigned int i = 0; i < rxn_val.size(); ++i)
        {
          if (x2_test[i] > upper) rxn_val[i] = -(upper-x2[i]-dt*diff2[i])/dt;
          else if (x2_test[i] < lower) rxn_val[i] = -(lower-x2[i]-dt*diff2[i])/dt;
        } */
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
