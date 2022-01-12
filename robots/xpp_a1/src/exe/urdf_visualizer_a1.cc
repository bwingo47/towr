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

#include <iostream>
#include <map>
#include <memory>
#include <string>

#include <ros/init.h>

#include <xpp_a1/inverse_kinematics_a1_4.h>
#include <xpp_msgs/topic_names.h>
#include <xpp_states/joints.h>
#include <xpp_states/endeffector_mappings.h>

#include <xpp_vis/cartesian_joint_converter.h>
#include <xpp_vis/urdf_visualizer.h>

using namespace xpp;
using namespace quad;

int main(int argc, char *argv[])
{
  ::ros::init(argc, argv, "urdf_visualizer_a1");

  const std::string joint_desired_a1 = "xpp/joint_a1_des";

/*
 * This class subscribes to a Cartesian robot state message and publishes
 * the joint state through inverse kinematics converted, joint state message. This
 * can then be used to visualize URDFs in RVIZ.
 * */
  auto a1_ik = std::make_shared<InverseKinematicsA1_4>();
  CartesianJointConverter inv_kin_converter(a1_ik,
					    xpp_msgs::robot_state_desired,
					    joint_desired_a1);

  // urdf joint names
  int n_ee = a1_ik->GetEECount();
  int n_j  = A1legJointCount;
  std::vector<UrdfVisualizer::URDFName> joint_names(n_ee*n_j);
  joint_names.at(n_j*LF + HAA) = "FL_hip_joint"; // Hip
  joint_names.at(n_j*LF + HFE) = "FL_thigh_joint"; // Upper leg
  joint_names.at(n_j*LF + KFE) = "FL_calf_joint"; // Lower leg
  joint_names.at(n_j*RF + HAA) = "FR_hip_joint";
  joint_names.at(n_j*RF + HFE) = "FR_thigh_joint";
  joint_names.at(n_j*RF + KFE) = "FR_calf_joint";
  joint_names.at(n_j*LH + HAA) = "RL_hip_joint";
  joint_names.at(n_j*LH + HFE) = "RL_thigh_joint";
  joint_names.at(n_j*LH + KFE) = "RL_calf_joint";
  joint_names.at(n_j*RH + HAA) = "RR_hip_joint";
  joint_names.at(n_j*RH + HFE) = "RR_thigh_joint";
  joint_names.at(n_j*RH + KFE) = "RR_calf_joint";

  // joint_names.at(n_j*LF + HAA) = "FL_hip_motor_2_chassis_joint"; // Hip
  // joint_names.at(n_j*LF + HFE) = "FL_upper_leg_2_hip_motor_joint"; // Upper leg
  // joint_names.at(n_j*LF + KFE) = "FL_lower_leg_2_upper_leg_joint"; // Lower leg
  // joint_names.at(n_j*RF + HAA) = "FR_hip_motor_2_chassis_joint";
  // joint_names.at(n_j*RF + HFE) = "FR_upper_leg_2_hip_motor_joint";
  // joint_names.at(n_j*RF + KFE) = "FR_lower_leg_2_upper_leg_joint";
  // joint_names.at(n_j*LH + HAA) = "RL_hip_motor_2_chassis_joint";
  // joint_names.at(n_j*LH + HFE) = "RL_upper_leg_2_hip_motor_joint";
  // joint_names.at(n_j*LH + KFE) = "RL_lower_leg_2_upper_leg_joint";
  // joint_names.at(n_j*RH + HAA) = "RR_hip_motor_2_chassis_joint";
  // joint_names.at(n_j*RH + HFE) = "RR_upper_leg_2_hip_motor_joint";
  // joint_names.at(n_j*RH + KFE) = "RR_lower_leg_2_upper_leg_joint";

/*
 * This class is responsible for converting an xpp_msgs::RobotStateJoint
 * robot message and publishing the corresponding RVIZ transforms.
 */
  std::string urdf = "a1_rviz_urdf_robot_description";
  UrdfVisualizer a1_desired(urdf, joint_names, "base", "world",
			     joint_desired_a1, "a1_des");

  ::ros::spin();

  return 1;
}