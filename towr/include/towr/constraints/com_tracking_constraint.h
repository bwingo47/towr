//
// Created by bruce on 6/10/21.
//

#ifndef TOWR_INCLUDE_TOWR_CONSTRAINTS_COM_TRACKING_CONSTRAINT_H_
#define TOWR_INCLUDE_TOWR_CONSTRAINTS_COM_TRACKING_CONSTRAINT_H_

#include <towr/variables/spline.h>
#include <towr/variables/spline_holder_extended.h>
#include <towr/variables/spline_holder.h>
#include <towr/variables/euler_converter.h>

#include <towr/models/kinematic_model.h>

#include <towr/constraints/time_discretization_constraint.h>

namespace towr {

class COMTrackingConstraint : public TimeDiscretizationConstraint {
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
  COMTrackingConstraint(std::vector<Vector3d> target, double T, double dt,
                        const SplineHolderExtended& spline_holder);
  COMTrackingConstraint(double t, double dt, const std::string &constraint_name);

  COMTrackingConstraint(std::vector<Vector3d> target, double T, double dt,
                        const SplineHolder& spline_holder);

  virtual ~COMTrackingConstraint() = default;


 private:
  NodeSpline::Ptr com_linear_;       ///< the linear position of the base/CoM.
  std::vector<Vector3d> target_;

  VecBound node_bounds_;


  // see TimeDiscretizationConstraint for documentation
  void UpdateConstraintAtInstance (double t, int k, VectorXd& g) const override;
  void UpdateBoundsAtInstance (double t, int k, VecBound&) const override;
  void UpdateJacobianAtInstance(double t, int k, std::string, Jacobian&) const override;

  int GetRow(int node, int dimension) const;
};

}

#endif //TOWR_INCLUDE_TOWR_CONSTRAINTS_COM_TRACKING_CONSTRAINT_H_
