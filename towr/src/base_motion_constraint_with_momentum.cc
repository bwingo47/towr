//
// Created by bruce on 3/29/21.
//

#include <towr/constraints/base_motion_constraint_with_momentum.h>
#include <towr/variables/variable_names.h>
#include <towr/variables/cartesian_dimensions.h>
#include <towr/variables/spline_holder.h>

namespace towr {


BaseMotionConstraintWithMomentum::BaseMotionConstraintWithMomentum (double T, double dt,
                                            const SplineHolderExtended& spline_holder)
    :TimeDiscretizationConstraint(T, dt, "baseMotion")
{
  com_linear_  = spline_holder.com_linear_;
  com_angular_ = spline_holder.com_angular_;

  double dev_rad = 0.05;
  node_bounds_.resize(k6D);
  node_bounds_.at(AX) = Bounds(-dev_rad, dev_rad);
  node_bounds_.at(AY) = Bounds(-dev_rad, dev_rad);
  node_bounds_.at(AZ) = ifopt::NoBound;//Bounds(-dev_rad, dev_rad);

  double z_init = com_linear_->GetPoint(0.0).p().z();
  node_bounds_.at(LX) = ifopt::NoBound;
  node_bounds_.at(LY) = ifopt::NoBound;//Bounds(-0.05, 0.05);
  node_bounds_.at(LZ) = Bounds(z_init-0.02, z_init+0.1); // allow to move dev_z cm up and down

  int n_constraints_per_node = node_bounds_.size();
  SetRows(GetNumberOfNodes()*n_constraints_per_node);
}

void
BaseMotionConstraintWithMomentum::UpdateConstraintAtInstance (double t, int k,
                                                  VectorXd& g) const
{
  g.middleRows(GetRow(k, LX), k3D) = com_linear_->GetPoint(t).p();
  g.middleRows(GetRow(k, AX), k3D) = com_angular_->GetPoint(t).p();
}

void
BaseMotionConstraintWithMomentum::UpdateBoundsAtInstance (double t, int k, VecBound& bounds) const
{
  for (int dim=0; dim<node_bounds_.size(); ++dim)
    bounds.at(GetRow(k,dim)) = node_bounds_.at(dim);
}

void
BaseMotionConstraintWithMomentum::UpdateJacobianAtInstance (double t, int k,
                                                std::string var_set,
                                                Jacobian& jac) const
{
  if (var_set == id::com_ang_nodes)
    jac.middleRows(GetRow(k,AX), k3D) = com_angular_->GetJacobianWrtNodes(t, kPos);

  if (var_set == id::com_lin_nodes)
    jac.middleRows(GetRow(k,LX), k3D) = com_linear_->GetJacobianWrtNodes(t, kPos);
}

int
BaseMotionConstraintWithMomentum::GetRow (int node, int dim) const
{
  return node*node_bounds_.size() + dim;
}

} /* namespace towr */