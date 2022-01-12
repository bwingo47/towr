//
// Created by bruce on 5/28/21.
//

#ifndef TOWR_INCLUDE_TOWR_COSTS_TIME_DISCRETIZATION_COST_H_
#define TOWR_INCLUDE_TOWR_COSTS_TIME_DISCRETIZATION_COST_H_

#include <string>
#include <vector>

#include <ifopt/cost_term.h>
#include <iostream>

namespace towr {

/**
 * @brief Add cost at discretized times along a trajectory.
 *
 * Often one want to add a cost at specific times t along
 * the trajectory. This class is responsible for building the overall
 * Jacobian from the individual Jacobians at each time instance.
 *
 * @ingroup Constraints
 */
class TimeDiscretizationCost : public ifopt::CostTerm {
 public:
  using VecTimes = std::vector<double>;
  using Bounds   = ifopt::Bounds;

  /**
   * @brief Constructs a constraint for ifopt.
   * @param T  The total duration of the trajectory.
   * @param dt The discretization interval at which each constraint is evaluated.
   * @param constraint_name  The name of the constraint.
   */
  TimeDiscretizationCost (double T, double dt, std::string cost_name);

  /**
   * @brief construct a constraint for ifopt.
   * @param dts  Time instances at which to evaluate the constraints.
   * @param name The name of the constraint.
   */
  TimeDiscretizationCost (const VecTimes& dts, std::string name);
  virtual ~TimeDiscretizationCost () = default;

//        double GetCost() const override;

  void FillJacobianBlock (std::string var_set, Jacobian&) const override;

 protected:
  int GetNumberOfNodes() const;
  VecTimes dts_; ///< times at which the constraint is evaluated.

 private:
  /**
   * @brief Sets Jacobian rows at a specific time t, corresponding to node k.
   * @param t  The time along the trajcetory to set the bounds.
   * @param k  The index of the time t, so t=k*dt
   * @param var_set The name of the ifopt variables currently being queried for.
   * @param[in/out] jac  The complete Jacobian, for which the corresponding
   *                     row and columns must be set.
   */
  virtual void UpdateJacobianAtInstance(double t, int k, std::string var_set,
                                        Jacobian& jac) const = 0;
};

} /* namespace towr */


#endif //TOWR_INCLUDE_TOWR_COSTS_TIME_DISCRETIZATION_COST_H_
