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

#include <iostream>

#include <towr/constraints/dynamic_constraint_with_momentum.h>

#include <towr/variables/variable_names.h>
#include <towr/variables/cartesian_dimensions.h>

namespace towr {

DynamicConstraintWithMomentum::DynamicConstraintWithMomentum (const SingleRigidBodyDynamicsWithMomentum::Ptr& m,
                                      double T, double dt,
                                      const SplineHolderExtended& spline_holder)
    :TimeDiscretizationConstraint(T, dt, "dynamic")
{
  model_ = m;

  // link with up-to-date spline variables
  com_linear_  = spline_holder.com_linear_;
  com_angular_ = EulerConverter(spline_holder.com_angular_);
  ellipsoid_axes_ = spline_holder.ellipsoid_axes_;
  ellipsoid_orient_ = EulerConverter(spline_holder.ellipsoid_orient_);
  angular_momentum_ = spline_holder.angular_momentum_;
  ee_forces_    = spline_holder.ee_force_;
  ee_motion_    = spline_holder.ee_motion_;

  SetRows(GetNumberOfNodes()*k9D);
}

int
DynamicConstraintWithMomentum::GetRow (int k, Dim9D dimension) const
{
  return k9D*k + dimension;
}

void
DynamicConstraintWithMomentum::UpdateConstraintAtInstance(double t, int k, VectorXd& g) const
{
  UpdateModel(t);
  g.segment(GetRow(k,aX), k9D) = model_->GetDynamicViolation();
}

void
DynamicConstraintWithMomentum::UpdateBoundsAtInstance(double t, int k, VecBound& bounds) const
{
  for (auto dim : AllDim9D)
    bounds.at(GetRow(k,dim)) = ifopt::BoundZero;
}

// TODO: Need to fix this method
void
DynamicConstraintWithMomentum::UpdateJacobianAtInstance(double t, int k, std::string var_set,
                                            Jacobian& jac) const
{
  UpdateModel(t);

  int n = jac.cols();
  Jacobian jac_model(k9D,n);

  // sensitivity of dynamic constraint w.r.t base variables.
  if (var_set == id::com_lin_nodes) {
    Jacobian jac_base_lin_pos = com_linear_->GetJacobianWrtNodes(t,kPos);
    Jacobian jac_base_lin_acc = com_linear_->GetJacobianWrtNodes(t,kAcc);

    jac_model = model_->GetJacobianWrtBaseLin(jac_base_lin_pos,
                                              jac_base_lin_acc);
  }

  if (var_set == id::com_ang_nodes) {
    jac_model = model_->GetJacobianWrtBaseAng(com_angular_, t);
  }

  if(var_set == id::ellipsoid_semi_axes_nodes) {
      Jacobian jac_ellipsoid_axes_pos = ellipsoid_axes_->GetJacobianWrtNodes(t, kPos);
      jac_model = model_->GetJacobianWrtEllipsoidAxes(jac_ellipsoid_axes_pos);
  }

  if(var_set == id::ellipsoid_orientation_nodes) {
      jac_model = model_->GetJacobianWrtEllipsoidOrient(ellipsoid_orient_, t);
  }

  if(var_set == id::angular_momentum_nodes) {
      Jacobian jac_angular_momentum_pos = angular_momentum_->GetJacobianWrtNodes(t, kPos);
      Jacobian jac_angular_momentum_vel = angular_momentum_->GetJacobianWrtNodes(t, kVel);
      jac_model = model_->GetJacobianWrtAM(jac_angular_momentum_pos, jac_angular_momentum_vel);
  }

  // sensitivity of dynamic constraint w.r.t. endeffector variables
  for (int ee=0; ee<model_->GetEECount(); ++ee) {
    if (var_set == id::EEForceNodes(ee)) {
      Jacobian jac_ee_force = ee_forces_.at(ee)->GetJacobianWrtNodes(t,kPos);
      jac_model = model_->GetJacobianWrtForce(jac_ee_force, ee);
    }

    if (var_set == id::EEMotionNodes(ee)) {
      Jacobian jac_ee_pos = ee_motion_.at(ee)->GetJacobianWrtNodes(t,kPos);
      jac_model = model_->GetJacobianWrtEEPos(jac_ee_pos, ee);
    }

    if (var_set == id::EESchedule(ee)) {
      Jacobian jac_f_dT = ee_forces_.at(ee)->GetJacobianOfPosWrtDurations(t);
      jac_model += model_->GetJacobianWrtForce(jac_f_dT, ee);

      Jacobian jac_x_dT = ee_motion_.at(ee)->GetJacobianOfPosWrtDurations(t);
      jac_model +=  model_->GetJacobianWrtEEPos(jac_x_dT, ee);
    }
  }

  jac.middleRows(GetRow(k,aX), k9D) = jac_model;
}

void
DynamicConstraintWithMomentum::UpdateModel (double t) const
{
  auto com = com_linear_->GetPoint(t);

  Eigen::Matrix3d w_R_b = com_angular_.GetRotationMatrixBaseToWorld(t);
  Eigen::Vector3d omega = com_angular_.GetAngularVelocityInWorld(t);
  Eigen::Vector3d omega_dot = com_angular_.GetAngularAccelerationInWorld(t);


  auto ellipsoid_axes_node = ellipsoid_axes_->GetPoint(t);
  // the 'e' vector
  auto ellipsoid_axes = ellipsoid_axes_node.p();
  // rotation matrix R(v) from base coord. to world coord.
  Eigen::Matrix3d ellipsoid_w_R_b = ellipsoid_orient_.GetRotationMatrixBaseToWorld(t);

  auto am_node = angular_momentum_->GetPoint(t);
  // angular momentum vector H
  auto am = am_node.p();
  // angular momentum rate vector H_dot
  auto am_dot = am_node.v();


  int n_ee = model_->GetEECount();
  std::vector<Eigen::Vector3d> ee_pos;
  std::vector<Eigen::Vector3d> ee_force;
  for (int ee=0; ee<n_ee; ++ee) {
//    std::cout << "ee=0" << ee_forces_.at(0)->GetPoint(t).p().transpose() << std::endl;
//    std::cout << "ee=1" << ee_forces_.at(1)->GetPoint(t).p().transpose() << std::endl;
//    std::cout << "ee=2" << ee_forces_.at(2)->GetPoint(t).p().transpose() << std::endl;
//    std::cout << "ee=3" << ee_forces_.at(3)->GetPoint(t).p().transpose() << std::endl;
    ee_force.push_back(ee_forces_.at(ee)->GetPoint(t).p());
    ee_pos.push_back(ee_motion_.at(ee)->GetPoint(t).p());
  }

  // TODO: DONE, check new implementation
  model_->SetCurrent(com.p(), com.a(), w_R_b, omega, omega_dot, ee_force, ee_pos,
                     ellipsoid_axes, ellipsoid_w_R_b, am, am_dot);
}

} /* namespace towr */
