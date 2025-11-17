// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attribute_loader.h>
#include <prismspf/core/variable_container.h>

#include <prismspf/config.h>
#include <prismspf/utilities/utilities.h>


PRISMS_PF_BEGIN_NAMESPACE
template <typename T, unsigned int dim> 
dealii::Tensor<2, voigt_tensor_size<dim>, T>
form_Cij(const dealii::Tensor<1, dim, T> &Cv1,
        const dealii::Tensor<1, dim, T> &Cv2,
        const dealii::Tensor<1, dim, T> &Cv3,
        const dealii::Tensor<1, dim, T> &Cv4,
        const dealii::Tensor<1, dim, T> &Cv5,
        const dealii::Tensor<1, dim, T> &Cv6,
        const dealii::Tensor<1, dim, T> &Cv7);

void
CustomAttributeLoader::load_variable_attributes()
{

  std::string c_vec_string = ", Cv1, Cv2, Cv3, Cv4, Cv5, Cv6, Cv7";
  set_variable_name(0, "u");
  set_variable_type(0, Vector);
  set_variable_equation_type(0, TimeIndependent);
  set_dependencies_gradient_term_rhs(0, std::string("grad(u)")+c_vec_string);
  set_dependencies_gradient_term_lhs(0, std::string("grad(change(u))")+c_vec_string);

  for (unsigned int i = 1; i <= 7; ++i)
    {
      set_variable_name(i, "Cv" + std::to_string(i));
      set_variable_type(i, Vector);
      set_variable_equation_type(i, Constant);
      set_dependencies_gradient_term_rhs(i, "");
      set_dependencies_gradient_term_lhs(i, "");
    }
  
  set_variable_name(8, "stress_diag"); 
  set_variable_type(8, Vector);
  set_variable_equation_type(8, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(8, std::string("grad(u)")+c_vec_string);
  set_is_postprocessed_field(8, true);

  set_variable_name(9, "stress_offdiag"); 
  set_variable_type(9, Vector);
  set_variable_equation_type(9, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(9, std::string("grad(u)")+c_vec_string);
  set_is_postprocessed_field(9, true);

  set_variable_name(10, "stress_eigs"); 
  set_variable_type(10, Vector);
  set_variable_equation_type(10, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(10, std::string("grad(u)")+c_vec_string);
  set_is_postprocessed_field(10, true);
  
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
  [[maybe_unused]] Types::Index                           index) const
{
  if (index == 0)
    {
      // Compute the RHS of the momentum balance equation
      dealii::Tensor<2, voigt_tensor_size<dim>, dealii::VectorizedArray<number>> 
      C = form_Cij<dealii::VectorizedArray<number>,dim>(
            variable_list.template get_value<ScalarGrad>(1),
            variable_list.template get_value<ScalarGrad>(2),
            variable_list.template get_value<ScalarGrad>(3),
            variable_list.template get_value<ScalarGrad>(4),
            variable_list.template get_value<ScalarGrad>(5),
            variable_list.template get_value<ScalarGrad>(6),
            variable_list.template get_value<ScalarGrad>(7));
      VectorGrad ux = variable_list.template get_symmetric_gradient<VectorGrad>(0);
      VectorGrad stress;
      compute_stress<dim, ScalarValue>(C, ux, stress);
      variable_list.set_gradient_term(0, -stress);
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
  if (index == 0)
    {
      dealii::Tensor<2, voigt_tensor_size<dim>, dealii::VectorizedArray<number>> 
      C = form_Cij<dealii::VectorizedArray<number>,dim>(
            variable_list.template get_value<ScalarGrad>(1),
            variable_list.template get_value<ScalarGrad>(2),
            variable_list.template get_value<ScalarGrad>(3),
            variable_list.template get_value<ScalarGrad>(4),
            variable_list.template get_value<ScalarGrad>(5),
            variable_list.template get_value<ScalarGrad>(6),
            variable_list.template get_value<ScalarGrad>(7));
      VectorGrad change_ux = variable_list.template get_symmetric_gradient<VectorGrad>(0,Change);
      VectorGrad stress;
      compute_stress<dim, ScalarValue>(C, change_ux, stress);
      variable_list.set_gradient_term(0, stress, Change);
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
    dealii::Tensor<2, voigt_tensor_size<dim>, dealii::VectorizedArray<number>> 
    C = form_Cij<dealii::VectorizedArray<number>,dim>(
          variable_list.template get_value<ScalarGrad>(1),
          variable_list.template get_value<ScalarGrad>(2),
          variable_list.template get_value<ScalarGrad>(3),
          variable_list.template get_value<ScalarGrad>(4),
          variable_list.template get_value<ScalarGrad>(5),
          variable_list.template get_value<ScalarGrad>(6),
          variable_list.template get_value<ScalarGrad>(7));
    VectorGrad ux = variable_list.template get_symmetric_gradient<VectorGrad>(0);
    VectorGrad stress;
    compute_stress<dim, ScalarValue>(C, ux, stress);
    ScalarGrad stress_diag;
    ScalarGrad stress_offdiag;
    for (unsigned int i = 0; i < dim; ++i)
      {
        stress_diag[i] = stress[i][i];
      }
    stress_offdiag[0] = stress[0][1];
    stress_offdiag[1] = stress[0][2];
    stress_offdiag[2] = stress[1][2];
    variable_list.set_value_term(8, stress_diag);
    variable_list.set_value_term(9, stress_offdiag);
    ScalarGrad stress_eigs;
    for (unsigned int k = 0; k < stress[0][0].size(); ++k)
      {
        dealii::SymmetricTensor<2, dim, number> stress_sym;
        for (unsigned int i = 0; i < dim; ++i)
          {
            for (unsigned int j = 0; j < dim; ++j)
              {
                stress_sym[i][j] = stress[i][j][k];
              }
          }
        std::array<number, dim> stress_eigvals = dealii::eigenvalues(stress_sym);
        for (unsigned int i = 0; i < dim; ++i)
          {
            stress_eigs[i][k] = stress_eigvals[i];
          }
      }
    variable_list.set_value_term(10, stress_eigs);
}

template <typename T, unsigned int dim> 
dealii::Tensor<2, voigt_tensor_size<dim>, T>
form_Cij(const dealii::Tensor<1, dim, T> &Cv1,
        const dealii::Tensor<1, dim, T> &Cv2,
        const dealii::Tensor<1, dim, T> &Cv3,
        const dealii::Tensor<1, dim, T> &Cv4,
        const dealii::Tensor<1, dim, T> &Cv5,
        const dealii::Tensor<1, dim, T> &Cv6,
        const dealii::Tensor<1, dim, T> &Cv7)
{
  if constexpr (dim != 3)
    {
      throw("This function is designed for 3-D elastic constants only.");
    }
  if constexpr (dim == 3)
    {
        dealii::Tensor<2, voigt_tensor_size<dim>, T> C;
        C[0][0] = Cv1[0];
        C[1][1] = Cv1[1];
        C[2][2] = Cv1[2];
        C[3][3] = Cv2[0];
        C[4][4] = Cv2[1];
        C[5][5] = Cv2[2];
        C[0][1] = Cv3[0];
        C[1][0] = Cv3[0];
        C[0][2] = Cv3[1];
        C[2][0] = Cv3[1];
        C[0][3] = Cv3[2];
        C[3][0] = Cv3[2];
        C[0][4] = Cv4[0];
        C[4][0] = Cv4[0];
        C[0][5] = Cv4[1];
        C[5][0] = Cv4[1];
        C[1][2] = Cv4[2];
        C[2][1] = Cv4[2];
        C[1][3] = Cv5[0];
        C[3][1] = Cv5[0];
        C[1][4] = Cv5[1];
        C[4][1] = Cv5[1];
        C[1][5] = Cv5[2];
        C[5][1] = Cv5[2];
        C[2][3] = Cv6[0];
        C[3][2] = Cv6[0];
        C[2][4] = Cv6[1];
        C[4][2] = Cv6[1];
        C[2][5] = Cv6[2];
        C[5][2] = Cv6[2];
        C[3][4] = Cv7[0];
        C[4][3] = Cv7[0];
        C[3][5] = Cv7[1];
        C[5][3] = Cv7[1];
        C[4][5] = Cv7[2];
        C[5][4] = Cv7[2];
        return C;
    }
}

#include "custom_pde.inst"

PRISMS_PF_END_NAMESPACE
