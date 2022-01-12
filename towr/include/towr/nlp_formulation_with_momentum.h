//
// Created by bruce on 3/18/21.
//

#ifndef TOWR_NLP_FORMULATION_WITH_MOMENTUM_H
#define TOWR_NLP_FORMULATION_WITH_MOMENTUM_H

#include <ifopt/variable_set.h>
#include <ifopt/constraint_set.h>
#include <ifopt/cost_term.h>

#include <towr/variables/spline_holder.h>
#include <towr/models/robot_model.h>
#include <towr/terrain/height_map.h>
#include <towr/parameters.h>
#include <towr/variables/spline_holder_extended.h>

#include <towr/costs/com_tracking_cost.h>

namespace towr {

/**
 * @defgroup Constraints
 * @brief Constraints of the trajectory optimization problem.
 *
 * These are the constraint sets that characterize legged locomotion.
 *
 * Folder: @ref include/towr/constraints
 */

/**
 * @defgroup Costs
 * @brief Costs of the trajectory optimization problem.
 *
 * These the the cost terms that prioritize certain solutions to the
 * legged locomotion problem.
 *
 * Folder: @ref include/towr/costs
 */

/**
 *
 * @brief A sample combination of variables, cost and constraints.
 *
 * This is _one_ example of how to combine the variables, constraints and costs
 * provided by this library. Additional variables or constraints can be added
 * to the NLP, or existing elements replaced to find a more powerful/general
 * formulation. This formulation was used to generate the motions described
 * in this paper: https://ieeexplore.ieee.org/document/8283570/
 */
class NlpFormulationMomentum {
public:
    using VariablePtrVec   = std::vector<ifopt::VariableSet::Ptr>;
    using ConstraintPtrVec  = std::vector<ifopt::ConstraintSet::Ptr>;
    using CostPtrVec       = std::vector<ifopt::CostTerm::Ptr>;
    using EEPos            = std::vector<Eigen::Vector3d>;
    using Vector3d         = Eigen::Vector3d;

    NlpFormulationMomentum ();
    virtual ~NlpFormulationMomentum () = default;

    /**
     * @brief The ifopt variable sets that will be optimized over.
     * @param[in/out] builds fully-constructed splines from the variables.
     */
    VariablePtrVec GetVariableSets(SplineHolderExtended& spline_holder_extended);

    /**
     * @brief The ifopt constraints that enforce feasible motions.
     * @param[in] uses the fully-constructed splines for initialization of constraints.
     */
    ConstraintPtrVec GetConstraints(const SplineHolderExtended& spline_holder_extended) const;

    /** @brief The ifopt costs to tune the motion. */
    ConstraintPtrVec GetCosts() const;

    ConstraintPtrVec GetTrackingCosts(const SplineHolderExtended& spline_holder_extended) const;

    ConstraintPtrVec GetNodeTrackingCosts() const;


//    BaseState initial_base_;
//    BaseState final_base_;
    CoMState initial_com_;
    CoMState final_com_;
//    Node initial_am_{3};
//    Node final_am_{3};
    DebugDummyState dummy_state1_;
    DebugDummyState dummy_state2_;
//    Node dummy_state3_;
    AngularMomentumState initial_am_;
    AngularMomentumState final_am_;
    EllipsoidState initial_ellipsoid_state_;
    EllipsoidState final_ellipsoid_state_;
    EEPos  initial_ee_W_;
    RobotModel model_;
    HeightMap::Ptr terrain_;
    Parameters params_;

private:
    // variables
    std::vector<NodesVariables::Ptr> MakeBaseVariables() const;
    std::vector<NodesVariables::Ptr> MakeCoMVariables() const;
    std::vector<NodesVariables::Ptr> MakeAngularMomentumVariables() const;
    std::vector<NodesVariables::Ptr> MakeEllipsoidVariables() const;
    std::vector<NodesVariablesPhaseBased::Ptr> MakeEndeffectorVariables() const;
    std::vector<NodesVariablesPhaseBased::Ptr> MakeForceVariables() const;
    std::vector<PhaseDurations::Ptr> MakeContactScheduleVariables() const;

    // constraints
    ConstraintPtrVec GetConstraint(Parameters::ConstraintName name,
                                  const SplineHolderExtended& splines) const;
    ConstraintPtrVec MakeDynamicConstraint(const SplineHolder& s) const;
    ConstraintPtrVec MakeDynamicConstraintWithMomentum(const SplineHolderExtended& s) const;
    ConstraintPtrVec MakeRangeOfMotionBoxConstraint(const SplineHolderExtended& s) const;
    ConstraintPtrVec MakeTotalTimeConstraint() const;
    ConstraintPtrVec MakeTerrainConstraint() const;
    ConstraintPtrVec MakeForceConstraint() const;
    ConstraintPtrVec MakeSwingConstraint() const;
    ConstraintPtrVec MakeBaseRangeOfMotionConstraint(const SplineHolderExtended& s) const;
    ConstraintPtrVec MakeBaseAccConstraint(const SplineHolderExtended& s) const;
    ConstraintPtrVec MakeCoMTrackingConstraint(const SplineHolderExtended& s) const;

    // costs
    CostPtrVec GetCost(const Parameters::CostName& id, double weight) const;
    CostPtrVec MakeForcesCost(double weight) const;
    CostPtrVec MakeEEMotionCost(double weight) const;

    CostPtrVec GetTrackingCost(const Parameters::CostName& id, double weight, std::vector<Eigen::Vector3d> target,
                               const SplineHolderExtended& s) const;
    CostPtrVec MakeEETrackingCost(std::vector<Eigen::Vector3d> target, double weight, const SplineHolderExtended& s) const;
    CostPtrVec MakeCOMTrackingCost(std::vector<Eigen::Vector3d> target, double weight, const SplineHolderExtended& s) const;
    CostPtrVec MakeAMTrackingCost(std::vector<Eigen::Vector3d> target, double weight, const SplineHolderExtended& s) const;
    CostPtrVec MakeLMTrackingCost(std::vector<Eigen::Vector3d> target, double weight, const SplineHolderExtended& s) const;
    CostPtrVec MakeMoITrackingCost(std::vector<Eigen::Vector3d> target, double weight, const SplineHolderExtended& s) const;

    CostPtrVec GetNodeTrackingCost(const Parameters::CostName& id, double weight, const std::vector<Eigen::Vector3d>& target) const;
    CostPtrVec MakeCOMNodeTrackingCost(const std::vector<Eigen::Vector3d>& target, double weight) const;

    // helper
    /**
     * @brief take a splineHolderExtended pointer and returns a pointer to a SplineHolder pointer,
     * removes the ellipsoid variables, and angular momentum variables.
     * @param s
     * @return
     */
    SplineHolder& SplineHolderExtendedToRegular(const SplineHolderExtended& s) const;
};

} /* namespace towr */


#endif //TOWR_NLP_FORMULATION_WITH_MOMENTUM_H
