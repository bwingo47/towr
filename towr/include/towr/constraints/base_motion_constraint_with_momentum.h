//
// Created by bruce on 3/29/21.
//

#ifndef TOWR_INCLUDE_TOWR_CONSTRAINTS_BASE_MOTION_CONSTRAINT_WITH_MOMENTUM_H_
#define TOWR_INCLUDE_TOWR_CONSTRAINTS_BASE_MOTION_CONSTRAINT_WITH_MOMENTUM_H_

#include <towr/variables/spline_holder.h>
#include <towr/variables/spline.h>
#include <towr/variables/spline_holder_extended.h>

#include "time_discretization_constraint.h"

namespace towr {

/**
 * @brief Keeps the 6D base motion in a specified range.
 *
 * In general this constraint should be avoided, since a similar affect can be
 * achieved with RangeOfMotionConstraint.
 *
 * @ingroup Constraints
 */
class BaseMotionConstraintWithMomentum : public TimeDiscretizationConstraint {
 public:
  /**
   * @brief Links the base variables and sets hardcoded bounds on the state.
   * @param T  The total time of the optimization horizon.
   * @param dt The discretization interval of the constraints.
   * @param spline_holder  Holds pointers to the base variables.
   */
  BaseMotionConstraintWithMomentum (double T, double dt, const SplineHolderExtended& spline_holder);
  virtual ~BaseMotionConstraintWithMomentum () = default;

  void UpdateConstraintAtInstance (double t, int k, VectorXd& g) const override;
  void UpdateBoundsAtInstance (double t, int k, VecBound&) const override;
  void UpdateJacobianAtInstance(double t, int k, std::string, Jacobian&) const override;

 private:
  NodeSpline::Ptr com_linear_;
  NodeSpline::Ptr com_angular_;

  VecBound node_bounds_;     ///< same bounds for each discretized node
  int GetRow (int node, int dim) const;
};

} /* namespace towr */

#endif //TOWR_INCLUDE_TOWR_CONSTRAINTS_BASE_MOTION_CONSTRAINT_WITH_MOMENTUM_H_
