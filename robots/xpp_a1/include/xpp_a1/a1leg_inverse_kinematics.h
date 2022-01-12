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

#ifndef XPP_VIS_A1LEG_INVERSE_KINEMATICS_H_
#define XPP_VIS_A1LEG_INVERSE_KINEMATICS_H_

#include <Eigen/Dense>
#include <ros/console.h>

namespace xpp {

enum A1JointID {HAA=0, HFE, KFE, A1legJointCount};

/**
 * @brief Converts a A1 foot position to joint angles.
 */
class A1legInverseKinematics {
public:
  using Vector3d = Eigen::Vector3d;
  enum KneeBend { Left, Right }; // used for symmetric leg, not for A1
  // mutable double q_aile_1, q_epaule_1, q_coude_1, q_aile_2, q_epaule_2, q_coude_2, q_aile_3, q_epaule_3, q_coude_3, q_aile_4, q_epaule_4, q_coude_4;
  // mutable int compteur ;
  /**
   * @brief Default c'tor initializing leg lengths with standard values.
   */
  A1legInverseKinematics () = default;
  virtual ~A1legInverseKinematics () = default;

  /**
   * @brief Returns the joint angles to reach a Cartesian foot position.
   * @param ee_pos_H  Foot position xyz expressed in the frame attached
   * at the hip-aa (H).
   */
  Vector3d GetJointAngles(const Vector3d& ee_pos_H, KneeBend bend) const;

  /**
   * @brief Restricts the joint angles to lie inside the feasible range
   * @param q[in/out]  Current joint angle that is adapted if it exceeds
   * the specified range.
   * @param joint  Which joint (HAA, HFE, KFE) this value represents.
   */
  void EnforceLimits(double& q, A1JointID joint) const;

private:
  // we base everything on the front-left leg for signs (LR)
  // Vector3d haa_to_hfe = Vector3d(0, 0.08505, 0); //distance of HAA to HFE
  // Vector3d hfe_to_kfe = Vector3d(0, 0, -0.2); //distance of HFE to KFE
  // Vector3d kfe_to_foot = Vector3d(0, 0, -0.2); // distance of KFE to FOOT

  Vector3d hfe_to_haa_y = Vector3d(0.0, 0.08505, 0); //distance of HFE to HAA in z direction
  double length_thigh = 0.2; // length of upper leg
  double length_shank = 0.2; // length of lower leg
};

} /* namespace xpp */

#endif /* XPP_VIS_A1LEG_INVERSE_KINEMATICS_H_ */