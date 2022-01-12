//
// Created by bruce on 5/31/21.
//

#ifndef TOWR_INCLUDE_TOWR_COSTS_NODETRACKINGCOST_H_
#define TOWR_INCLUDE_TOWR_COSTS_NODETRACKINGCOST_H_

#include <memory>
#include <string>

#include <Eigen/Dense>

#include <ifopt/cost_term.h>

#include <towr/variables/nodes_variables.h>

namespace towr {

class NodeTrackingCost : public ifopt::CostTerm{
 public:

  NodeTrackingCost(const std::string& nodes_id, Dx deriv, int dim, const std::vector<Eigen::Vector3d>& target, double weight);
  virtual ~NodeTrackingCost() = default;

  void InitVariableDependedQuantities(const VariablesPtr& x) override;

  double GetCost() const override;

 private:

  std::shared_ptr<NodesVariables> nodes_;

  std::vector<Eigen::Vector3d> target_;

  std::string node_id_;
  Dx deriv_;
  int dim_;
  double weight_;

  void FillJacobianBlock(std::string var_set, Jacobian&) const override;

};


}

#endif //TOWR_INCLUDE_TOWR_COSTS_NODETRACKINGCOST_H_
