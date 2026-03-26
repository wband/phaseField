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
  set_variable_type(0, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(0, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(0, "c, n, mu");
  set_solve_block(0, 1);

  set_variable_name(1, "u");
  set_variable_type(1, FieldInfo::TensorRank::Vector);
  set_variable_equation_type(1, TimeIndependent);
  set_dependencies_value_term_rhs(1, "");
  set_dependencies_gradient_term_rhs(1, dependency_string);
  set_dependencies_value_term_lhs(1, "");
  set_dependencies_gradient_term_lhs(1, "c, n, grad(change(u))");
  set_solve_block(1, 0);

  set_variable_name(2, "mu");
  set_variable_type(2, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(2, Auxiliary);
  set_dependencies_value_term_rhs(2, "c, n, grad(n), grad(c), grad(u)");
  set_dependencies_gradient_term_rhs(2, "n, grad(n), grad(c)");
  set_solve_block(2, 0);

  set_variable_name(3, "n");
  set_variable_type(3, FieldInfo::TensorRank::Scalar);
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
      ScalarValue c   = variable_list.template get_value<ScalarValue>(0);
      ScalarGrad  cx  = variable_list.template get_gradient<ScalarGrad>(0);
      ScalarValue n  = variable_list.template get_value<ScalarValue>(1);
      VectorGrad  ux  = variable_list.template get_symmetric_gradient<VectorGrad>(2);
      
      VectorGrad e_el = ux - sfts * n;
      VectorGrad stress, stress_c_effect;
      dealii::Tensor<2, voigt_tensor_size<dim>, ScalarValue> CIJ_combined =
         CIJ_FP * (1.0 - c) + CIJ_LFP * c;
      compute_stress<dim, ScalarValue>(CIJ_combined, strain, stress);
      compute_stress<dim, ScalarValue>(CIJ_LFP-CIJ_FP, strain, stress_c_effect);
      
      ScalarValue elastic_driving_force = 0.0;
      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              elastic_driving_force += 0.5 * stress[i][j] * sfts_const1[i][j];
              elastic_driving_force += 0.5 * stress_c_effect[i][j] * sfts_const1[i][j];
            }
         }
      }

  // Compute strain
  VectorGrad strain = ux - (sfts_const1 * h1V + sfts_const2 * h2V + sfts_const3 * h3V);

  // Compute stress
  VectorGrad stress;

  dealii::Tensor<2, voigt_tensor_size<dim>, ScalarValue> CIJ_combined =
    CIJ_Mg * (1.0 - sum_hV) + CIJ_Beta * sum_hV;
  compute_stress<dim, ScalarValue>(CIJ_combined, strain, stress);

  // Compute one of the stress terms in the order parameter chemical potential,
  // nDependentMisfitACp = C*(E-E0)*(E0_p*Hn)
  ScalarValue nDependentMisfitAC1 = 0.0;
  ScalarValue nDependentMisfitAC2 = 0.0;
  ScalarValue nDependentMisfitAC3 = 0.0;
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          nDependentMisfitAC1 += stress[i][j] * sfts_const1[i][j];
          nDependentMisfitAC2 += stress[i][j] * sfts_const2[i][j];
          nDependentMisfitAC3 += stress[i][j] * sfts_const3[i][j];
        }
    }

  nDependentMisfitAC1 *= -hn1V;
  nDependentMisfitAC2 *= -hn2V;
  nDependentMisfitAC3 *= -hn3V;

  // Compute the other stress term in the order parameter chemical potential,
  // heterMechACp = 0.5*Hn*(C_beta-C_alpha)*(E-E0)*(E-E0)
  ScalarValue heterMechAC1 = 0.0;
  ScalarValue heterMechAC2 = 0.0;
  ScalarValue heterMechAC3 = 0.0;
  if (n_dependent_stiffness == true)
    {
      VectorGrad stress_2;
      compute_stress<dim, ScalarValue>(CIJ_Beta - CIJ_Mg, strain, stress_2);
      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              heterMechAC1 += stress_2[i][j] * strain[i][j];
            }
        }
      // Aside from HnpV, heterMechAC1, heterMechAC2, and heterMechAC3 are equal
      heterMechAC2 = 0.5 * hn2V * heterMechAC1;
      heterMechAC3 = 0.5 * hn3V * heterMechAC1;
      heterMechAC1 = 0.5 * hn1V * heterMechAC1;
    }

  // compute K*nx
  ScalarGrad Knx1;
  for (unsigned int a = 0; a < dim; a++)
    {
      Knx1[a] = 0.0;
      for (unsigned int b = 0; b < dim; b++)
        {
          Kcx[a] += Kcxij[a][b] * cx[b];
        }
    }

  ScalarValue sum_hV = h1V + h2V + h3V;

  // The terms in the govering equations
  ScalarValue eq_c       = c;
  ScalarGrad  eqx_c_temp = cx * ((1.0 - sum_hV) * faccV + sum_hV * fbccV) +
                          n1x * ((fbcV - facV) * hn1V) + n2x * ((fbcV - facV) * hn2V) +
                          n3x * ((fbcV - facV) * hn3V);
  ScalarGrad  eqx_c = -get_timestep() * McV * eqx_c_temp;
  ScalarValue eq_n1 = n1 - get_timestep() * Mn1V *
                             ((fbV - faV) * hn1V + nDependentMisfitAC1 + heterMechAC1);
  ScalarValue eq_n2 = n2 - get_timestep() * Mn2V *
                             ((fbV - faV) * hn2V + nDependentMisfitAC2 + heterMechAC2);
  ScalarValue eq_n3 = n3 - get_timestep() * Mn3V *
                             ((fbV - faV) * hn3V + nDependentMisfitAC3 + heterMechAC3);
  ScalarGrad eqx_n1 = -get_timestep() * Mn1V * Knx1;
  ScalarGrad eqx_n2 = -get_timestep() * Mn2V * Knx2;
  ScalarGrad eqx_n3 = -get_timestep() * Mn3V * Knx3;

  variable_list.set_value_term(0, eq_c);
  variable_list.set_gradient_term(0, eqx_c);
  variable_list.set_value_term(1, eq_n1);
  variable_list.set_gradient_term(1, eqx_n1);
  variable_list.set_value_term(2, eq_n2);
  variable_list.set_gradient_term(2, eqx_n2);
  variable_list.set_value_term(3, eq_n3);
  variable_list.set_gradient_term(3, eqx_n3);
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
  if (index == 4)
    {
      ScalarValue n1 = variable_list.template get_value<ScalarValue>(1);
      ScalarValue n2 = variable_list.template get_value<ScalarValue>(2);
      ScalarValue n3 = variable_list.template get_value<ScalarValue>(3);
      VectorGrad  ux = variable_list.template get_symmetric_gradient<VectorGrad>(4);

      // Interpolation functions
      ScalarValue h1V = compute_hV(n1);
      ScalarValue h2V = compute_hV(n2);
      ScalarValue h3V = compute_hV(n3);

      // Compute strain
      VectorGrad strain =
        ux - (sfts_const1 * h1V + sfts_const2 * h2V + sfts_const3 * h3V);

      // Compute stress
      VectorGrad stress;
      if (n_dependent_stiffness == true)
        {
          ScalarValue                                            sum_hV = h1V + h2V + h3V;
          dealii::Tensor<2, voigt_tensor_size<dim>, ScalarValue> CIJ_combined =
            CIJ_Mg * (1.0 - sum_hV) + CIJ_Beta * sum_hV;
          compute_stress<dim, ScalarValue>(CIJ_combined, strain, stress);
        }
      else
        {
          compute_stress<dim, ScalarValue>(CIJ_Mg, strain, stress);
        }

      variable_list.set_gradient_term(4, -stress);
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
{
  if (index == 4)
    {
      ScalarValue n1 = variable_list.template get_value<ScalarValue>(1);
      ScalarValue n2 = variable_list.template get_value<ScalarValue>(2);
      ScalarValue n3 = variable_list.template get_value<ScalarValue>(3);
      VectorGrad  change_ux =
        variable_list.template get_symmetric_gradient<VectorGrad>(4, Change);

      // Interpolation functions
      ScalarValue h1V = compute_hV(n1);
      ScalarValue h2V = compute_hV(n2);
      ScalarValue h3V = compute_hV(n3);

      // Compute strain
      VectorGrad strain = change_ux;

      // Compute stress
      VectorGrad stress;
      if (n_dependent_stiffness == true)
        {
          ScalarValue                                            sum_hV = h1V + h2V + h3V;
          dealii::Tensor<2, voigt_tensor_size<dim>, ScalarValue> CIJ_combined =
            CIJ_Mg * (1.0 - sum_hV) + CIJ_Beta * sum_hV;
          compute_stress<dim, ScalarValue>(CIJ_combined, strain, stress);
        }
      else
        {
          compute_stress<dim, ScalarValue>(CIJ_Mg, strain, stress);
        }

      variable_list.set_gradient_term(4, stress, Change);
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
  ScalarValue c   = variable_list.template get_value<ScalarValue>(0);
  ScalarValue n1  = variable_list.template get_value<ScalarValue>(1);
  ScalarGrad  n1x = variable_list.template get_gradient<ScalarGrad>(1);
  ScalarValue n2  = variable_list.template get_value<ScalarValue>(2);
  ScalarGrad  n2x = variable_list.template get_gradient<ScalarGrad>(2);
  ScalarValue n3  = variable_list.template get_value<ScalarValue>(3);
  ScalarGrad  n3x = variable_list.template get_gradient<ScalarGrad>(3);
  VectorGrad  ux  = variable_list.template get_symmetric_gradient<VectorGrad>(4);

  // Free energy expressions and interpolation functions
  ScalarValue faV    = A0 + A1 * c + A2 * c * c + A3 * c * c * c + A4 * c * c * c * c;
  ScalarValue fbV    = B2 * c * c + B1 * c + B0;
  ScalarValue h1V    = compute_hV(n1);
  ScalarValue h2V    = compute_hV(n2);
  ScalarValue h3V    = compute_hV(n3);
  ScalarValue sum_hV = h1V + h2V + h3V;

  ScalarValue f_chem = (1.0 - sum_hV) * faV + sum_hV * fbV;

  ScalarValue f_grad = 0.0;
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          f_grad += 0.5 * Kn1[i][j] * n1x[i] * n1x[j] +
                    0.5 * Kn2[i][j] * n2x[i] * n2x[j] + 0.5 * Kn3[i][j] * n3x[i] * n3x[j];
        }
    }

  // Compute strain
  VectorGrad strain = ux - (sfts_const1 * h1V + sfts_const2 * h2V + sfts_const3 * h3V);

  // Compute stress
  VectorGrad stress;
  if (n_dependent_stiffness == true)
    {
      dealii::Tensor<2, voigt_tensor_size<dim>, ScalarValue> CIJ_combined =
        CIJ_Mg * (1.0 - sum_hV) + CIJ_Beta * sum_hV;
      compute_stress<dim, ScalarValue>(CIJ_combined, strain, stress);
    }
  else
    {
      compute_stress<dim, ScalarValue>(CIJ_Mg, strain, stress);
    }

  ScalarValue f_el = 0.0;
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          f_el += 0.5 * stress[i][j] * strain[i][j];
        }
    }

  ScalarValue f_tot = f_chem + f_grad + f_el;

  variable_list.set_value_term(5, f_tot);
}

#include "custom_pde.inst"

PRISMS_PF_END_NAMESPACE
