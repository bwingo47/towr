//
// Created by bruce on 6/9/21.
//

#include "towr/costs/ee_tracking_cost.h"
#include "towr/variables/variable_names.h"

namespace towr {

EETrackingCost::EETrackingCost (std::vector<Vector3d> target, double weight,
                                double T, double dt,
                                const EE& ee,
                                const SplineHolderExtended& spline_holder)
    :TimeDiscretizationCost(T, dt, "ee-tracking-cost-" + std::to_string(ee))
{
  // No need to specify the number of rows because there is only one row for the cost
  ee_motion_    = spline_holder.ee_motion_.at(ee);
  ee_ = ee;
  target_ = target;
  weight_ = weight;
}

EETrackingCost::EETrackingCost (std::vector<Vector3d> target, double weight,
                                double T, double dt,
                                const EE& ee,
                                const SplineHolder& spline_holder)
    :TimeDiscretizationCost(T, dt, "ee-tracking-cost-" + std::to_string(ee))
{
  // No need to specify the number of rows because there is only one row for the cost
  ee_motion_    = spline_holder.ee_motion_.at(ee);
  ee_ = ee;
  target_ = target;
  weight_ = weight;
}

int
EETrackingCost::GetRow (int node, int dim) const
{
  return node*k3D + dim;
}

double
EETrackingCost::GetCost() const {
  double cost;
//  assert(target_.size() == GetNumberOfNodes());

  int k = 0;
  for (double t : dts_) {
    Vector3d pos_ee_W = ee_motion_->GetPoint(t).p();
    Vector3d tar = target_.at(k);
    Vector3d ee_error = pos_ee_W - tar;
    cost += weight_*ee_error.transpose()*ee_error;
    k++;
  }

  return cost;
}

void
EETrackingCost::UpdateJacobianAtInstance (double t, int k,
                                          std::string var_set,
                                          Jacobian& jac) const
{
  Vector3d pos_ee_W = ee_motion_->GetPoint(t).p();
//        int row_start = GetRow(k,X);

  // dc/dt = 2*(EE(t) - ref)*(dEE(t)/dEEMotionNodes)
  if (var_set == id::EEMotionNodes(ee_)) {
    jac.middleRows(0, 1) += 2*weight_*(pos_ee_W - target_.at(k)).transpose().sparseView()
        *ee_motion_->GetJacobianWrtNodes(t,kPos);
  }

}

} /* namespace admm */

