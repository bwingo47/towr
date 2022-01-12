//
// Created by bruce on 6/10/21.
//

#include <towr/constraints/com_tracking_constraint.h>
#include <towr/variables/cartesian_dimensions.h>
#include <towr/variables/variable_names.h>

namespace towr {

COMTrackingConstraint::COMTrackingConstraint(std::vector<Vector3d> target, double T, double dt,
                                             const SplineHolderExtended& spline_holder)
    : TimeDiscretizationConstraint(T, dt, "comTracking")
{
  com_linear_ = spline_holder.com_linear_;
  target_ = target;

  double dev_x = 0.05;
  double dev_y = 0.05;
  double dev_z = 0.3;
  node_bounds_.resize(k3D);
  node_bounds_.at(X) = Bounds(-dev_x, dev_x);
  node_bounds_.at(Y) = Bounds(-dev_y, dev_y);
  node_bounds_.at(Z) = Bounds(-dev_z, dev_z);

  int n_constraints_per_node = node_bounds_.size();
  SetRows(GetNumberOfNodes()*n_constraints_per_node);
}

COMTrackingConstraint::COMTrackingConstraint(std::vector<Vector3d> target, double T, double dt,
                                             const SplineHolder& spline_holder)
    : TimeDiscretizationConstraint(T, dt, "comTracking")
{
  com_linear_ = spline_holder.base_linear_;
  target_ = target;
  double dev_x = 0.05;
  double dev_y = 0.05;
  double dev_z = 0.3;
  node_bounds_.resize(k3D);
  node_bounds_.at(X) = Bounds(-dev_x, dev_x);
  node_bounds_.at(Y) = Bounds(-dev_y, dev_y);
  node_bounds_.at(Z) = Bounds(-dev_z, dev_z);

  int n_constraints_per_node = node_bounds_.size();
  SetRows(GetNumberOfNodes()*n_constraints_per_node);
}

int
COMTrackingConstraint::GetRow(int node, int dimension) const
{
  return node*node_bounds_.size() + dimension;
}

void
COMTrackingConstraint::UpdateConstraintAtInstance(double t, int k, VectorXd& g) const
{
  Vector3d com_linear_curr = com_linear_->GetPoint(t).p();
  Vector3d com_target = target_[k];
  Vector3d com_err = com_linear_curr - com_target;

  // Hardmard product
//  g.middleRows(GetRow(k, X), k3D) = com_err.cwiseProduct(com_err);
  g.middleRows(GetRow(k, X), k3D) = com_err;
}

void
COMTrackingConstraint::UpdateBoundsAtInstance(double t, int k, VecBound& bounds) const
{
  for (int dim=0; dim<node_bounds_.size(); ++dim) {
    bounds.at(GetRow(k,dim)) = node_bounds_.at(dim);
  }
}

void
COMTrackingConstraint::UpdateJacobianAtInstance(double t, int k,
                                                std::string var_set,
                                                Jacobian& jac) const
{
  // derivative of Hardmard product:
  // replace the Hadamard product of a vector by a matrix product with the diagonal matrix formed from that vector
  // e.g. f = x O y ==> f = X*y, where O is Hardmard product, and * is matrix product.
  // X = Diag(x)

  Vector3d com_linear_curr = com_linear_->GetPoint(t).p();
  Vector3d com_target = target_.at(k);
  Vector3d com_err = com_linear_curr - com_target;

  Eigen::Matrix3d com_err_matrix = com_err.asDiagonal();

  int row_start = GetRow(k,X);

  if (var_set == id::base_lin_nodes) {
//    jac.middleRows(row_start, k3D) = com_err_matrix.sparseView()*com_linear_->GetJacobianWrtNodes(t, kPos);
    jac.middleRows(row_start, k3D) = com_linear_->GetJacobianWrtNodes(t, kPos);
  }

}



}