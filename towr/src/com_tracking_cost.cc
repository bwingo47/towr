//
// Created by bruce on 5/28/21.
//


#include "towr/costs/com_tracking_cost.h"
#include <towr/variables/variable_names.h>

namespace towr {

COMTrackingCost::COMTrackingCost (std::vector<Vector3d> target, double weight,
                                  double T, double dt,
                                  const SplineHolderExtended& spline_holder)
    :TimeDiscretizationCost(T, dt, "com-tracking-cost")
{
  // No need to specify the number of rows because there is only one row for the cost
  com_linear_  = spline_holder.com_linear_;
  target_ = target;
  weight_ = weight;
}

COMTrackingCost::COMTrackingCost (std::vector<Vector3d> target, double weight,
                                  double T, double dt,
                                  const SplineHolder& spline_holder)
    :TimeDiscretizationCost(T, dt, "com-tracking-cost")
{
  // No need to specify the number of rows because there is only one row for the cost
  com_linear_  = spline_holder.base_linear_;
  target_ = target;
  weight_ = weight;
}

int
COMTrackingCost::GetRow (int node, int dim) const
{
  return node*k3D + dim;
}

double
COMTrackingCost::GetCost() const {
  double cost;
//  assert(target_.size() == GetNumberOfNodes());
  int k = 0;
  for (double t : dts_) {
    Vector3d base_pos_t = com_linear_->GetPoint(t).p();
    Vector3d tar = target_.at(k);
    Vector3d com_error = base_pos_t - tar;
    cost += weight_*com_error.transpose()*com_error;
    k++;

//            GetVariables()->GetComponent(id::com_lin_nodes)->GetValues();
  }

  return cost;
}

void
COMTrackingCost::UpdateJacobianAtInstance (double t, int k,
                                           std::string var_set,
                                           Jacobian& jac) const
{
  Vector3d base_pos_t = com_linear_->GetPoint(t).p();
//        int row_start = GetRow(k,X);

  // dc/dt = 2*(EE(t) - ref)*(dEE(t)/dNodes)
  if ((var_set == id::com_lin_nodes) || (var_set == id::base_lin_nodes)) {
    jac.middleRows(0, 1) += 2*weight_*(base_pos_t - target_.at(k)).transpose().sparseView()
        *com_linear_->GetJacobianWrtNodes(t, kPos);

//            Eigen::SparseMatrix<double, Eigen::RowMajor> jac_debug = com_linear_->GetJacobianWrtNodes(t, kPos);
//    Eigen::SparseMatrix<double, Eigen::RowMajor> jac_debug = 2*weight_*(base_pos_t - target_.at(k)).transpose().sparseView()*com_linear_->GetJacobianWrtNodes(t, kPos);
//    Eigen::MatrixXd jac_debug_dense;
//    jac_debug_dense = Eigen::MatrixXd(jac_debug);
//    std::cout << " " << std::endl;
  }

}

} /* namespace towr */


