// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/initial_conditions.h>
#include <prismspf/core/nonuniform_dirichlet.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/utilities/utilities.h>

#include <prismspf/config.h>

#include <algorithm>
#include <numbers>
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
  using std::max;
  using std::min;
  using std::sin;
  using std::numbers::pi;
  auto   domain_size = get_user_inputs().get_spatial_discretization().get_size();
  
  if (index == 0)
    {
      double y = point[1];
      double x = point[0];
      double ys = (0.5 * (domain_size[1]*0.75 - y
       + std::sin(2.0*3.14*(2.845 * x + 1.0)/domain_size[0]) * 2.0 
       + std::sin(2.0*3.14*(7.123 * x      )/domain_size[0]) * 1.0));
    //  double ys = (0.5 * (domain_size[1]*0.75 - y));
      double flat = 0.5 * (1.0 + sin(pi * max(-0.5, min(0.5, std::sqrt(2.0) * ys /l_int))));
      unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

 /*     // Declaring distribution (Type std::uniform_real_distribution<double>)
      std::uniform_real_distribution<double> dist(
        std::uniform_real_distribution<double>(0.0, 1.0));

      // Declaring random number generator (Type std::mt19937_64)
      auto rng = std::mt19937_64(seed);
      flat = flat + 2.0 * 0.001 * (dist(rng) - 0.5); */

      scalar_value = max(min(flat,1.0-1e-4),1e-4);
    }
  else if (index == 1)
    {
      scalar_value = 0.2;
    }
  else if (index == 2)
    {
      scalar_value = 0.05;
    }
  else 
    {
      scalar_value = 0.0;
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
