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
  [[maybe_unused]] number                   &scalar_value,
  [[maybe_unused]] number                   &vector_component_value) const
{
  const number ylength = this->get_user_inputs().get_spatial_discretization().get_size()[1];
  number posx   = point[0];
  number posy   = point[1];
  number channel_width = ylength/4;
  scalar_value = 0.0;
  if(index == 3)
    {
      if (posx < channel_depth)
        {
          scalar_value = 0.5 * (1.0 + tanh((abs(posy - ylength/2.0) - channel_width/2)/ell));
        }
      else 
        {
          number dist = sqrt((posx - channel_depth)*(posx - channel_depth) + (posy - ylength/2.0)*(posy - ylength/2.0));
          scalar_value = 0.5 * (1.0 + tanh((dist - channel_width/2.0)/ell));
        }
      scalar_value = std::max<number>(scalar_value, 1e-6);
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

#include "custom_pde.inst"

PRISMS_PF_END_NAMESPACE
