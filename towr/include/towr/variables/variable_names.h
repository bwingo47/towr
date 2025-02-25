/******************************************************************************
Copyright (c) 2018, Alexander W. Winkler. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
******************************************************************************/

#ifndef TOWR_VARIABLES_VARIABLE_NAMES_H_
#define TOWR_VARIABLES_VARIABLE_NAMES_H_

#include <string>

namespace towr {
/**
 * @brief Identifiers (names) used for variables in the optimization problem.
 *
 * @ingroup Variables
 */
namespace id {

static const std::string base_lin_nodes    = "base-lin";
static const std::string base_ang_nodes    = "base-ang";
static const std::string com_lin_nodes    = "com-lin";
static const std::string angular_momentum_nodes    = "angular-momentum";
static const std::string com_ang_nodes    = "com-ang";
static const std::string ellipsoid_semi_axes_nodes   = "ellip-axes";
static const std::string ellipsoid_orientation_nodes   = "ellip-orientation";

static const std::string ee_motion_nodes   = "ee-motion_";
static const std::string ee_force_nodes    = "ee-force_";
static const std::string angular_momentum_phase_based_nodes    = "angular-momentum-phase-based";
static const std::string contact_schedule  = "ee-schedule";


static std::string EEMotionNodes(uint ee)
{
  return  ee_motion_nodes + std::to_string(ee);
}

static std::string EEForceNodes(uint ee)
{
  return  ee_force_nodes + std::to_string(ee);
}

static std::string EESchedule(uint ee)
{
  return  contact_schedule + std::to_string(ee);
}

} // namespace id
} // namespace towr



#endif /* TOWR_VARIABLES_VARIABLE_NAMES_H_ */
