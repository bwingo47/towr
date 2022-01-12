//
// Created by bruce on 5/31/21.
//

#include <towr/costs/node_tracking_cost.h>
#include <cmath>

namespace towr {

NodeTrackingCost::NodeTrackingCost(const std::string &nodes_id,
                                   towr::Dx deriv,
                                   int dim,
                                   const std::vector<Eigen::Vector3d> &target,
                                   double weight)
                                   : CostTerm(nodes_id
                                   +"-dx_"+std::to_string(deriv)
                                   +"-dim_"+std::to_string(dim))
{
  node_id_ = nodes_id;
  weight_ = weight;
  deriv_ = deriv;
  dim_ = dim;
  target_ = target;
}

void
NodeTrackingCost::InitVariableDependedQuantities(const ifopt::ConstraintSet::VariablesPtr &x)
{
  nodes_ = x->GetComponent<NodesVariables>(node_id_);
}

double
NodeTrackingCost::GetCost() const
{
  double cost;
  int num_nodes = (nodes_->GetNodes()).size();
  int num_target = target_.size();
  int increment = floor(num_target/num_nodes);
  int idx = 0;

  for (auto n : nodes_->GetNodes()) {
    double val = n.at(deriv_)(dim_);
    double target = target_.at(idx*increment)(dim_);
    cost += weight_*std::pow(val-target,2);
    idx++;
  }
  return cost;
}

void
NodeTrackingCost::FillJacobianBlock(std::string var_set, ifopt::Component::Jacobian& jac) const
{
  int num_nodes = (nodes_->GetNodes()).size();
  int num_target = target_.size();
  int increment = floor(num_target/num_nodes);

  if (var_set == node_id_) {
    for (int i=0; i<nodes_->GetRows(); ++i)
      for (auto nvi : nodes_->GetNodeValuesInfo(i))
        if (nvi.deriv_==deriv_ && nvi.dim_==dim_) {
          double val = nodes_->GetNodes().at(nvi.id_).at(deriv_)(dim_);
//          double target = target_.at(i*increment)(dim_);
//          double target = 2.0;
          jac.coeffRef(0, i) += weight_*2.0*(val-target_.at(i*increment)(dim_));
        }
  }

}

}

