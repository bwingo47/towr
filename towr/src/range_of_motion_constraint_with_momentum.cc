//
// Created by bruce on 3/29/21.
//

#include <towr/constraints/range_of_motion_constraint_with_momentum.h>
#include <towr/variables/variable_names.h>

namespace towr {

RangeOfMotionConstraintWithMomentum::RangeOfMotionConstraintWithMomentum (const KinematicModel::Ptr& model,
                                                                          double T, double dt,
                                                                          const EE& ee,
                                                                          const SplineHolderExtended& spline_holder)
                                                                         :TimeDiscretizationConstraint(T, dt, "rangeofmotion-" + std::to_string(ee))
{
  com_linear_  = spline_holder.com_linear_;
  com_angular_ = EulerConverter(spline_holder.com_angular_);
  ee_motion_    = spline_holder.ee_motion_.at(ee);

  max_deviation_from_nominal_ = model->GetMaximumDeviationFromNominal();
  nominal_ee_pos_B_           = model->GetNominalStanceInBase().at(ee);
  ee_ = ee;

  SetRows(GetNumberOfNodes()*k3D);
}

int
RangeOfMotionConstraintWithMomentum::GetRow (int node, int dim) const
{
  return node*k3D + dim;
}

void
RangeOfMotionConstraintWithMomentum::UpdateConstraintAtInstance (double t, int k, VectorXd& g) const
{
  Vector3d base_W  = com_linear_->GetPoint(t).p();
  Vector3d pos_ee_W = ee_motion_->GetPoint(t).p();
  EulerConverter::MatrixSXd b_R_w = com_angular_.GetRotationMatrixBaseToWorld(t).transpose();

  Vector3d vector_base_to_ee_W = pos_ee_W - base_W;
  Vector3d vector_base_to_ee_B = b_R_w*(vector_base_to_ee_W);

  g.middleRows(GetRow(k, X), k3D) = vector_base_to_ee_B;
}

void
RangeOfMotionConstraintWithMomentum::UpdateBoundsAtInstance (double t, int k, VecBound& bounds) const
{
  for (int dim=0; dim<k3D; ++dim) {
    ifopt::Bounds b;
    b += nominal_ee_pos_B_(dim);
    b.upper_ += max_deviation_from_nominal_(dim);
    b.lower_ -= max_deviation_from_nominal_(dim);
    bounds.at(GetRow(k,dim)) = b;
  }
}

void
RangeOfMotionConstraintWithMomentum::UpdateJacobianAtInstance (double t, int k,
                                                   std::string var_set,
                                                   Jacobian& jac) const
{
  EulerConverter::MatrixSXd b_R_w = com_angular_.GetRotationMatrixBaseToWorld(t).transpose();
  int row_start = GetRow(k,X);

  if (var_set == id::com_lin_nodes) {
    jac.middleRows(row_start, k3D) = -1*b_R_w*com_linear_->GetJacobianWrtNodes(t, kPos);
  }

  if (var_set == id::com_ang_nodes) {
    Vector3d base_W   = com_linear_->GetPoint(t).p();
    Vector3d ee_pos_W = ee_motion_->GetPoint(t).p();
    Vector3d r_W = ee_pos_W - base_W;
    jac.middleRows(row_start, k3D) = com_angular_.DerivOfRotVecMult(t,r_W, true);
  }

  if (var_set == id::EEMotionNodes(ee_)) {
    jac.middleRows(row_start, k3D) = b_R_w*ee_motion_->GetJacobianWrtNodes(t,kPos);
  }

  if (var_set == id::EESchedule(ee_)) {
    jac.middleRows(row_start, k3D) = b_R_w*ee_motion_->GetJacobianOfPosWrtDurations(t);
  }
}

} /* namespace xpp */