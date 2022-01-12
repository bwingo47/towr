//
// Created by bruce on 6/9/21.
//

#ifndef TOWR_INCLUDE_TOWR_COSTS_EE_TRACKING_COST_H_
#define TOWR_INCLUDE_TOWR_COSTS_EE_TRACKING_COST_H_

#include <towr/variables/spline.h>
#include <towr/variables/spline_holder_extended.h>
#include <towr/variables/spline_holder.h>
#include <towr/variables/euler_converter.h>

#include <towr/models/kinematic_model.h>

#include <towr/costs/time_discretization_cost.h>

namespace towr {

/** @brief Constrains an endeffector to lie in a box around the nominal stance.
  *
  * These constraints are necessary to avoid configurations
  * that are outside the kinematic reach of the robot. The constraint
  * is defined by Cartesian estimates of the reachability of each endeffector.
  *
  * This constraint calculates the position of of the contact expressed in the
  * current CoM frame and constrains it to lie in a box around the nominal/
  * natural contact position for that leg.
  *
  * @ingroup Constraints
  */
class EETrackingCost : public TimeDiscretizationCost {
 public:
  using EE = uint;
  using Vector3d = Eigen::Vector3d;

  /**
   * @brief Constructs a constraint instance.
   * @param robot_model   The kinematic restrictions of the robot.
   * @param T   The total duration of the optimization.
   * @param dt  the discretization intervall at which to enforce constraints.
   * @param ee            The endeffector for which to constrain the range.
   * @param spline_holder Pointer to the current variables.
   */
  EETrackingCost(std::vector<Vector3d> target, double weight, double T, double dt,
                 const EE& ee,
                 const SplineHolderExtended& spline_holder);

  EETrackingCost(std::vector<Vector3d> target, double weight, double T, double dt,
                 const EE& ee,
                 const SplineHolder& spline_holder);

  virtual ~EETrackingCost() = default;
  double GetCost () const override;

 private:
  NodeSpline::Ptr ee_motion_;       ///< the linear position of the endeffectors.
  EE ee_;
  std::vector<Vector3d> target_;
  double weight_;

  // see TimeDiscretizationConstraint for documentation
  void UpdateJacobianAtInstance(double t, int k, std::string, Jacobian&) const override;

  int GetRow(int node, int dimension) const;
};

} /* namespace towr */

#endif //TOWR_INCLUDE_TOWR_COSTS_EE_TRACKING_COST_H_
