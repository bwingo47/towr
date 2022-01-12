/******************************************************************************
Copyright (c) 2017, Alexander W. Winkler. All rights reserved.
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

#include <xpp_a1/a1leg_inverse_kinematics.h>

#include <cmath>
#include <map>

#include <xpp_states/cartesian_declarations.h>


namespace xpp {


A1legInverseKinematics::Vector3d
A1legInverseKinematics::GetJointAngles (const Vector3d& ee_pos_B, KneeBend bend) const
{
  double q_HAA, q_HFE, q_KFE; //hip, shoulder, knee
  Eigen::Vector3d xr;
  Eigen::Matrix3d R;

  // translate to the local coordinate of the attachment of the leg
  // and flip coordinate signs such that all computations can be done
  // for the front-left leg
  xr = ee_pos_B;

  // compute the HAA angle
  if (bend==Left)
    q_HAA = atan2(xr[Y],-xr[Z])-atan2(0.08505,sqrt(pow(xr[Z],2)+pow(xr[Y],2)-pow(0.08505,2)));
  else
    q_HAA = -atan2(xr[Y],-xr[Z])+atan2(0.08505,sqrt(pow(xr[Z],2)+pow(xr[Y],2)-pow(0.08505,2)));
  
  // rotate into the HFE coordinate system (rot around X)
  R << 1.0, 0.0, 0.0, 0.0, cos(q_HAA), -sin(q_HAA), 0.0, sin(q_HAA), cos(q_HAA);

  xr = (R * xr).eval();

  // translate into the HFE coordinate system (along Z axis)
  xr += hfe_to_haa_y;  //distance of HFE to HAA in z direction

  // compute square of length from HFE to foot
  double tmp1 = pow(xr[X],2)+pow(xr[Z],2);


  // compute temporary angles (with reachability check)
  double lu = length_thigh;  // length of upper leg
  double ll = length_shank;  // length of lower leg
  double alpha = atan2(-xr[Z],xr[X]) - 0.5*M_PI;  //  flip and rotate to match HyQ joint definition


  double some_random_value_for_beta = (pow(lu,2)+tmp1-pow(ll,2))/(2.*lu*sqrt(tmp1)); // this must be between -1 and 1
  if (some_random_value_for_beta > 1) {
    some_random_value_for_beta = 1;
  }
  if (some_random_value_for_beta < -1) {
    some_random_value_for_beta = -1;
  }
  double beta = acos(some_random_value_for_beta);

  // compute Hip FE angle
  q_HFE = alpha + beta;


  double some_random_value_for_gamma = (pow(ll,2)+pow(lu,2)-tmp1)/(2.*ll*lu);
  // law of cosines give the knee angle
  if (some_random_value_for_gamma > 1) {
    some_random_value_for_gamma = 1;
  }
  if (some_random_value_for_gamma < -1) {
    some_random_value_for_gamma = -1;
  }
  double gamma  = acos(some_random_value_for_gamma);

  q_KFE = gamma - M_PI;

  // knee bend
  EnforceLimits(q_HAA, HAA);
  EnforceLimits(q_HFE, HFE);
  EnforceLimits(q_KFE, KFE);

  // if (bend==Forward)
  //   return Vector3d(q_HAA_bf, q_HFE_bf, q_KFE_bf);
  // else // backward
  return Vector3d(q_HAA, q_HFE, q_KFE);

}

void
A1legInverseKinematics::EnforceLimits (double& val, A1JointID joint) const
{
  // joint limit for A1
  // const static double haa_min = -0.872664625997;
  // const static double haa_max =  1.0471975512;

  // const static double hfe_min = -0.523598775598;
  // const static double hfe_max =  3.92699081699;

  // const static double kfe_min = -2.77507351067;
  // const static double kfe_max = -0.610865238198;

  // totally exaggerated joint angle limits
  const static double haa_min = -M_PI/2;
  const static double haa_max =  M_PI/4;

  const static double hfe_min = -M_PI/4;
  const static double hfe_max =  M_PI/4;

  const static double kfe_min = -M_PI/2;
  const static double kfe_max =  0;

  // reduced joint angles for optimization
  static const std::map<A1JointID, double> max_range {
    {HAA, haa_max},
    {HFE, hfe_max},
    {KFE, kfe_max}
  };

  // reduced joint angles for optimization
  static const std::map<A1JointID, double> min_range {
    {HAA, haa_min},
    {HFE, hfe_min},
    {KFE, kfe_min}
  };

  double max = max_range.at(joint);
  val = val>max? max : val;

  double min = min_range.at(joint);
  val = val<min? min : val;
}

} /* namespace xpp */