// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/initial_conditions.h>
#include <prismspf/core/nonuniform_dirichlet.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

#include <cmath>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::set_initial_condition(
  [[maybe_unused]] const unsigned int       &index,
  [[maybe_unused]] const unsigned int       &component,
  [[maybe_unused]] const dealii::Point<dim> &point,
  [[maybe_unused]] double                   &scalar_value,
  [[maybe_unused]] double                   &vector_component_value) const
{
     const double SYSTEM_WIDTH = this->get_user_inputs().get_spatial_discretization().get_size()[0];
      const double SYSTEM_HEIGHT = this->get_user_inputs().get_spatial_discretization().get_size()[1];
      double posx   = point[0];
      double posy   = point[1];
      double channel_width = SYSTEM_WIDTH / 4;
      double depth    = 0.7 * SYSTEM_HEIGHT;
      double vacancy_end = 0.3 * SYSTEM_HEIGHT;
      double vacancy_start = 0.75 * SYSTEM_HEIGHT;

      // double dist;
      // dist = 0.0;
      // for (unsigned int dir = 0; dir < dim; dir++){
      //     dist += (point[dir]-center[dir])*(point[dir]-center[dir]);
      // }
      // dist = std::sqrt(dist); 

      double channel_Left   = 0.5 * (1 + tanh((abs(posx - 1 * SYSTEM_WIDTH) - 5 * channel_width / 2) / ell));
      double channel_Right  = 0.5 * (1 - tanh((abs(posx - 1 * SYSTEM_WIDTH) - 3 * channel_width / 2) / ell));
      double channel_Bottom = 0.5 * (1 - tanh((posy - depth) / ell));
      double vacancy_Top = 0.5 * (1 + tanh((posy - vacancy_start) / ell));
      double vacancy_Bottom = 0.5 * (1 - tanh((posy - vacancy_end) / ell));
 
  if (index == 2){
      // scalar_value = mu_bar/kwell;
    // if (posx < 0.07 && posx > 0.05 && posy > 0.02 && posy < 0.04)
    //   scalar_value = n_init;
    // else
    //   scalar_value = 0.0;

    scalar_value = n_init*(1.0-(channel_Left+vacancy_Top+channel_Right+vacancy_Bottom))*std::max(posy - vacancy_end,0.0)/(depth-vacancy_end);
    if (scalar_value > n_init)
      scalar_value = n_init;
    if (scalar_value < 0.0)
      scalar_value = 0.0;
  }

  if (index ==1){
    vector_component_value = 0.0;
  }
  if (index == 0){

    // scalar_value = 1.0-(slabPsi + 1.0-(channel_Left+channel_Right+channel_diag));
    scalar_value = channel_Left+channel_Right+channel_Bottom;
    if (scalar_value > 1.0)
      scalar_value = 1.0;
    if (scalar_value <= 0.0)
      scalar_value = 1e-6;
  }
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::set_nonuniform_dirichlet(
  [[maybe_unused]] const unsigned int       &index,
  [[maybe_unused]] const unsigned int       &boundary_id,
  [[maybe_unused]] const unsigned int       &component,
  [[maybe_unused]] const dealii::Point<dim> &point,
  [[maybe_unused]] number                   &scalar_value,
  [[maybe_unused]] number                   &vector_component_value) const
{}

INSTANTIATE_TRI_TEMPLATE(CustomPDE)

PRISMS_PF_END_NAMESPACE
