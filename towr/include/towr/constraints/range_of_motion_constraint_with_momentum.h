//
// Created by bruce on 3/29/21.
//

#ifndef TOWR_INCLUDE_TOWR_CONSTRAINTS_RANGE_OF_MOTION_CONSTRAINT_WITH_MOMENTUM_H_
#define TOWR_INCLUDE_TOWR_CONSTRAINTS_RANGE_OF_MOTION_CONSTRAINT_WITH_MOMENTUM_H_

#include <towr/variables/spline.h>
#include <towr/variables/spline_holder.h>
#include <towr/variables/euler_converter.h>

#include <towr/models/kinematic_model.h>
#include <towr/variables/spline_holder_extended.h>

#include "time_discretization_constraint.h"

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
class RangeOfMotionConstraintWithMomentum : public TimeDiscretizationConstraint {
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
  RangeOfMotionConstraintWithMomentum(const KinematicModel::Ptr& robot_model,
                          double T, double dt,
                          const EE& ee,
                          const SplineHolderExtended& spline_holder);
  virtual ~RangeOfMotionConstraintWithMomentum() = default;

 private:
  NodeSpline::Ptr com_linear_;     ///< the linear position of the base.
  EulerConverter com_angular_; ///< the orientation of the base.
  NodeSpline::Ptr ee_motion_;       ///< the linear position of the endeffectors.

  Eigen::Vector3d max_deviation_from_nominal_;
  Eigen::Vector3d nominal_ee_pos_B_;
  EE ee_;

  // see TimeDiscretizationConstraint for documentation
  void UpdateConstraintAtInstance (double t, int k, VectorXd& g) const override;
  void UpdateBoundsAtInstance (double t, int k, VecBound&) const override;
  void UpdateJacobianAtInstance(double t, int k, std::string, Jacobian&) const override;

  int GetRow(int node, int dimension) const;
};

} /* namespace towr */

#endif //TOWR_INCLUDE_TOWR_CONSTRAINTS_RANGE_OF_MOTION_CONSTRAINT_WITH_MOMENTUM_H_
