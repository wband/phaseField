// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/pde_operator_base.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class CustomPDE : public PDEOperatorBase<dim, degree, number>
{
public:
  using ScalarValue = dealii::VectorizedArray<number>;
  using ScalarGrad  = dealii::Tensor<1, dim, ScalarValue>;
  using ScalarHess  = dealii::Tensor<2, dim, ScalarValue>;
  using VectorValue = dealii::Tensor<1, dim, ScalarValue>;
  using VectorGrad  = dealii::Tensor<2, dim, ScalarValue>;
  using VectorHess  = dealii::Tensor<3, dim, ScalarValue>;
  using PDEOperatorBase<dim, degree, number>::get_user_inputs;
  using PDEOperatorBase<dim, degree, number>::get_pf_tools;

  /**
   * @brief Constructor.
   */
  CustomPDE(const UserInputParameters<dim> &_user_inputs, PhaseFieldTools<dim> &_pf_tools)
    : PDEOperatorBase<dim, degree, number>(_user_inputs, _pf_tools)
    , m_well(get_user_inputs().user_constants.get_double("m_well"))
    , kappa(get_user_inputs().user_constants.get_double("kappa"))
  {}

private:

  void
  post_solve_block([[maybe_unused]] SolveContext<dim, degree, number> &solve_context,
                   [[maybe_unused]] unsigned int                       solver_id) override
  {
    
    if (solver_id == 2) // explicit n
      {
        number increment = solve_context.get_simulation_timer().get_increment();
        if (increment > -1)
          {
            const auto &n_test1_field =
                solve_context.get_solution_indexer().get_solution_vector(1);
            number n_test1 =
              n_test1_field * solve_context.get_invm_manager().get_jxw(TensorRank::Scalar);
            
            const auto &argument_field =
              solve_context.get_solution_indexer().get_solution_vector(2);
            dealii::LinearAlgebra::distributed::Vector<number> n_test2_field;
            n_test2_field.reinit(argument_field);
            for (unsigned int i = 0; i < n_test2_field.locally_owned_size(); ++i)
              {
                n_test2_field.local_element(i) =
                    std::sinh(argument_field.local_element(i));
              }
            number n_test2 = n_test2_field * solve_context.get_invm_manager().get_jxw(TensorRank::Scalar);
            Logger::instance() << LogFormatter::info("Test1 value: ") << std::setprecision(12) << n_test1 << "\n";
            Logger::instance() << LogFormatter::info("Test2 value: ") << std::setprecision(12) << n_test2 << "\n";
          }
      }
  }



  void
  set_initial_condition([[maybe_unused]] const unsigned int       &index,
                        [[maybe_unused]] const unsigned int       &component,
                        [[maybe_unused]] const dealii::Point<dim> &point,
                        [[maybe_unused]] number                   &scalar_value,
                        [[maybe_unused]] number &vector_component_value) const override
  {
    const dealii::Tensor<1, dim> &mesh_size =
      get_user_inputs().spatial_discretization.rectangular_mesh.size;
    if (index == 0)
      {
        scalar_value = 1.0  +  std::sin(point[0] * 3.14159 * 2.0);
      }
    if (index == 3)
      {
        scalar_value = std::sinh(std::abs(std::cos(point[0] * 3.14159 * 2.0)) * 3.14159 * 2.0);
      }
  }

  void
  compute_rhs(FieldContainer<dim, degree, number> &variable_list,
              const SimulationTimer               &sim_timer,
              unsigned int                         solve_block_id) const override
  {
    if (solve_block_id == 1) // explicit n
      {
        ScalarValue n_val  = variable_list.template get_value<Scalar, OldOne>(0);
        variable_list.set_value_term(0, n_val + 0.1);
      }
    else if (solve_block_id == 2)
      {
        ScalarValue n_val  = variable_list.template get_value<Scalar, Current>(0);
        ScalarGrad n_grad  = variable_list.template get_gradient<Scalar, Current>(0);
        ScalarValue argument = n_grad.norm();
        ScalarValue n_test1_val = std::sinh(argument);
        variable_list.set_value_term(1, n_test1_val);
        variable_list.set_value_term(2, argument);
      }
    else if (solve_block_id == 3)
      {
        ScalarValue argument = variable_list.template get_value<Scalar, Current>(2);
        variable_list.set_value_term(4, std::sinh(argument));
      }
  }

  number m_well;
  number kappa;
};

PRISMS_PF_END_NAMESPACE
