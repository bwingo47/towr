//
// Created by bruce on 5/28/21.
//

#ifndef TOWR_INCLUDE_TOWR_COSTS_COM_TRACKING_COST_H_
#define TOWR_INCLUDE_TOWR_COSTS_COM_TRACKING_COST_H_

#include <towr/variables/spline.h>
#include <towr/variables/spline_holder_extended.h>
#include <towr/variables/spline_holder.h>
#include <towr/variables/euler_converter.h>

#include <towr/models/kinematic_model.h>

#include <towr/costs/time_discretization_cost.h>

using namespace towr;

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
class COMTrackingCost : public TimeDiscretizationCost {
 public:
  using Vector3d = Eigen::Vector3d;

  /**
   * @brief Constructs a constraint instance.
   * @param robot_model   The kinematic restrictions of the robot.
   * @param T   The total duration of the optimization.
   * @param dt  the discretization intervall at which to enforce constraints.
   * @param ee            The endeffector for which to constrain the range.
   * @param spline_holder Pointer to the current variables.
   */
  COMTrackingCost(std::vector<Vector3d> target, double weight, double T, double dt,
                  const SplineHolderExtended& spline_holder);

  COMTrackingCost(std::vector<Vector3d> target, double weight, double T, double dt,
                  const SplineHolder& spline_holder);

  virtual ~COMTrackingCost() = default;
  double GetCost () const override;

 private:
  NodeSpline::Ptr com_linear_;       ///< the linear position of the base/CoM.
  std::vector<Vector3d> target_;
  double weight_;

  // see TimeDiscretizationConstraint for documentation
  void UpdateJacobianAtInstance(double t, int k, std::string, Jacobian&) const override;

  int GetRow(int node, int dimension) const;
};

} /* namespace towr */

#endif //TOWR_INCLUDE_TOWR_COSTS_COM_TRACKING_COST_H_
