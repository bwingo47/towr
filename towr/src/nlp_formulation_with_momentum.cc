//
// Created by bruce on 3/18/21.
//

#include <towr/nlp_formulation_with_momentum.h>

#include <towr/nlp_formulation.h>

#include <towr/variables/variable_names.h>
#include <towr/variables/phase_durations.h>

#include <towr/constraints/base_motion_constraint_with_momentum.h>
#include <towr/constraints/dynamic_constraint.h>
#include <towr/constraints/dynamic_constraint_with_momentum.h>
#include <towr/constraints/force_constraint.h>
#include <towr/constraints/range_of_motion_constraint_with_momentum.h>
#include <towr/constraints/swing_constraint.h>
#include <towr/constraints/terrain_constraint.h>
#include <towr/constraints/total_duration_constraint.h>
#include <towr/constraints/spline_acc_constraint.h>
#include <towr/constraints/com_tracking_constraint.h>

#include <towr/costs/node_cost.h>
#include <towr/costs/node_tracking_cost.h>
#include <towr/variables/nodes_variables_all.h>

#include <iostream>
#include <towr/variables/spline_holder_extended.h>

namespace towr {

    NlpFormulationMomentum::NlpFormulationMomentum() {
        using namespace std;
        cout << "\n";
        cout << "************************************************************\n";
        cout << " TOWR - Trajectory Optimization for Walking Robots (v1.4)\n";
        cout << "                \u00a9 Alexander W. Winkler\n";
        cout << "           https://github.com/ethz-adrl/towr\n";
        cout << "************************************************************";
        cout << "\n\n";
//        int a = -1;
    }

    NlpFormulationMomentum::VariablePtrVec
    NlpFormulationMomentum::GetVariableSets(SplineHolderExtended &spline_holder) {
        VariablePtrVec vars;

        auto com_motion = MakeCoMVariables();
        vars.insert(vars.end(), com_motion.begin(), com_motion.end());

        // make ellipsoid variables
        auto ellipsoid_motion = MakeEllipsoidVariables();
        vars.insert(vars.end(), ellipsoid_motion.begin(), ellipsoid_motion.end());

        // make angular momentum variables
        auto am_trajectory = MakeAngularMomentumVariables();
        vars.insert(vars.end(), am_trajectory.begin(), am_trajectory.end());

        auto ee_motion = MakeEndeffectorVariables();
        vars.insert(vars.end(), ee_motion.begin(), ee_motion.end());

        auto ee_force = MakeForceVariables();
        vars.insert(vars.end(), ee_force.begin(), ee_force.end());

        auto contact_schedule = MakeContactScheduleVariables();
        // can also just be fixed timings that aren't optimized over, but still added
        // to spline_holder.
        if (params_.IsOptimizeTimings()) {
            vars.insert(vars.end(), contact_schedule.begin(), contact_schedule.end());
        }


        // stores these readily constructed spline
        spline_holder = SplineHolderExtended(com_motion.at(0), // linear
                                     com_motion.at(1), // angular
                                     params_.GetCoMPolyDurations(),
                                     ellipsoid_motion.at(0), // principal semi-axes
                                     ellipsoid_motion.at(1), // orientation
                                     params_.GetEllipsoidPolyDurations(),
                                     am_trajectory.at(0), // angular momentum
                                     params_.GetAMPolyDurations(),
                                     ee_motion,
                                     ee_force,
                                     contact_schedule,
                                     params_.IsOptimizeTimings());
        return vars;
    }

    std::vector<NodesVariables::Ptr>
    NlpFormulationMomentum::MakeBaseVariables() const {
//        std::vector<NodesVariables::Ptr> vars;
//
//        int n_nodes = params_.GetBasePolyDurations().size() + 1;
//
//        auto spline_lin = std::make_shared<NodesVariablesAll>(n_nodes, k3D, id::base_lin_nodes);
//
//        double x = final_base_.lin.p().x();
//        double y = final_base_.lin.p().y();
//        double z = terrain_->GetHeight(x, y) - model_.kinematic_model_->GetNominalStanceInBase().front().z();
//        Vector3d final_pos(x, y, z);
//
//        spline_lin->SetByLinearInterpolation(initial_base_.lin.p(), final_pos, params_.GetTotalTime());
//        spline_lin->AddStartBound(kPos, {X, Y, Z}, initial_base_.lin.p());
//        spline_lin->AddStartBound(kVel, {X, Y, Z}, initial_base_.lin.v());
//        spline_lin->AddFinalBound(kPos, params_.bounds_final_lin_pos_, final_base_.lin.p());
//        spline_lin->AddFinalBound(kVel, params_.bounds_final_lin_vel_, final_base_.lin.v());
//        vars.push_back(spline_lin);
//
//        auto spline_ang = std::make_shared<NodesVariablesAll>(n_nodes, k3D, id::base_ang_nodes);
//        spline_ang->SetByLinearInterpolation(initial_base_.ang.p(), final_base_.ang.p(), params_.GetTotalTime());
//        spline_ang->AddStartBound(kPos, {X, Y, Z}, initial_base_.ang.p());
//        spline_ang->AddStartBound(kVel, {X, Y, Z}, initial_base_.ang.v());
//        spline_ang->AddFinalBound(kPos, params_.bounds_final_ang_pos_, final_base_.ang.p());
//        spline_ang->AddFinalBound(kVel, params_.bounds_final_ang_vel_, final_base_.ang.v());
//        vars.push_back(spline_ang);
//
//        return vars;
    }

    std::vector<NodesVariables::Ptr>
    NlpFormulationMomentum::MakeCoMVariables() const {
        std::vector<NodesVariables::Ptr> vars;

        int n_nodes = params_.GetCoMPolyDurations().size() + 1;

        auto spline_lin = std::make_shared<NodesVariablesAll>(n_nodes, k3D, id::com_lin_nodes);

        spline_lin->SetByLinearInterpolation(initial_com_.lin.p(), final_com_.lin.p(), params_.GetTotalTime());
        spline_lin->AddStartBound(kPos, {X, Y, Z}, initial_com_.lin.p());
        spline_lin->AddStartBound(kVel, {X, Y, Z}, initial_com_.lin.v());
        spline_lin->AddFinalBound(kPos, {X, Y, Z}, final_com_.lin.p());
        spline_lin->AddFinalBound(kVel, {X, Y, Z}, final_com_.lin.v());
        vars.push_back(spline_lin);

        auto spline_ang = std::make_shared<NodesVariablesAll>(n_nodes, k3D, id::com_ang_nodes);
        spline_ang->SetByLinearInterpolation(initial_com_.ang.p(), final_com_.ang.p(), params_.GetTotalTime());
        spline_ang->AddStartBound(kPos, {X, Y, Z}, initial_com_.ang.p());
        spline_ang->AddStartBound(kVel, {X, Y, Z}, initial_com_.ang.v());
        spline_ang->AddFinalBound(kPos, params_.bounds_final_ang_pos_, final_com_.ang.p());
        spline_ang->AddFinalBound(kVel, params_.bounds_final_ang_vel_, final_com_.ang.v());
        vars.push_back(spline_ang);

        return vars;
    }

    std::vector<NodesVariables::Ptr>
    NlpFormulationMomentum::MakeEllipsoidVariables() const {
        std::vector<NodesVariables::Ptr> vars;

        int n_nodes = params_.GetEllipsoidPolyDurations().size() + 1;

        auto spline_ellipsoid_axes = std::make_shared<NodesVariablesAll>(n_nodes, k3D,
                                                                         id::ellipsoid_semi_axes_nodes);
        spline_ellipsoid_axes->SetByLinearInterpolation(initial_ellipsoid_state_.axis.p(),
                                                        final_ellipsoid_state_.axis.p(), params_.GetTotalTime());
        spline_ellipsoid_axes->AddStartBound(kPos, {X, Y, Z}, initial_ellipsoid_state_.axis.p());
        spline_ellipsoid_axes->AddStartBound(kVel, {X, Y, Z}, initial_ellipsoid_state_.axis.v());
        spline_ellipsoid_axes->AddFinalBound(kPos, {X, Y, Z}, final_ellipsoid_state_.axis.p());
        spline_ellipsoid_axes->AddFinalBound(kVel, {X, Y, Z}, final_ellipsoid_state_.axis.v());
        vars.push_back(spline_ellipsoid_axes);

        auto spline_ellipsoid_orients = std::make_shared<NodesVariablesAll>(n_nodes, k3D,
                                                                            id::ellipsoid_orientation_nodes);
        spline_ellipsoid_orients->SetByLinearInterpolation(initial_ellipsoid_state_.orient.p(),
                                                           final_ellipsoid_state_.orient.p(),
                                                           params_.GetTotalTime());
        spline_ellipsoid_orients->AddStartBound(kPos, {X, Y, Z}, initial_ellipsoid_state_.orient.p());
        spline_ellipsoid_orients->AddStartBound(kVel, {X, Y, Z}, initial_ellipsoid_state_.orient.v());
        spline_ellipsoid_orients->AddFinalBound(kPos, {X, Y, Z}, final_ellipsoid_state_.orient.p());
        spline_ellipsoid_orients->AddFinalBound(kVel, {X, Y, Z}, final_ellipsoid_state_.orient.v());
        vars.push_back(spline_ellipsoid_orients);

        return vars;
    }

    std::vector<NodesVariables::Ptr>
    NlpFormulationMomentum::MakeAngularMomentumVariables() const {
//        std::vector<NodesVariables::Ptr> vars;
//
//        int n_nodes = params_.GetAMPolyDurations().size() + 1;
//
//        auto spline_am = std::make_shared<NodesVariablesAll>(n_nodes, k3D, id::angular_momentum_nodes);
//
//        spline_am->SetByLinearInterpolation(initial_am_.p(), final_am_.p(), params_.GetTotalTime());
//        spline_am->AddStartBound(kPos, {X, Y, Z}, initial_am_.p());
//        spline_am->AddStartBound(kVel, {X, Y, Z}, initial_am_.v());
//        spline_am->AddFinalBound(kPos, {X, Y, Z}, final_am_.p());
//        spline_am->AddFinalBound(kVel, {X, Y, Z}, final_am_.v());
//        vars.push_back(spline_am);
//
//        return vars;
      std::vector<NodesVariables::Ptr> vars;

      int n_nodes = params_.GetAMPolyDurations().size() + 1;

      auto spline_am = std::make_shared<NodesVariablesAll>(n_nodes, k3D, id::angular_momentum_nodes);

      spline_am->SetByLinearInterpolation(initial_am_.values.p(), final_am_.values.p(), params_.GetTotalTime());
      spline_am->AddStartBound(kPos, {X, Y, Z}, initial_am_.values.p());
      spline_am->AddStartBound(kVel, {X, Y, Z}, initial_am_.values.v());
      spline_am->AddFinalBound(kPos, {X, Y, Z}, final_am_.values.p());
      spline_am->AddFinalBound(kVel, {X, Y, Z}, final_am_.values.v());
      vars.push_back(spline_am);

      return vars;
    }


    std::vector<NodesVariablesPhaseBased::Ptr>
    NlpFormulationMomentum::MakeEndeffectorVariables() const {
        std::vector<NodesVariablesPhaseBased::Ptr> vars;

        // Endeffector Motions
        double T = params_.GetTotalTime();
        for (int ee = 0; ee < params_.GetEECount(); ee++) {
            auto nodes = std::make_shared<NodesVariablesEEMotion>(
                    params_.GetPhaseCount(ee),
                    params_.ee_in_contact_at_start_.at(ee),
                    id::EEMotionNodes(ee),
                    params_.ee_polynomials_per_swing_phase_);

            // initialize towards final footholds
            double yaw = final_com_.ang.p().z();
            Eigen::Vector3d euler(0.0, 0.0, yaw);
            Eigen::Matrix3d w_R_b = EulerConverter::GetRotationMatrixBaseToWorld(euler);
            Vector3d final_ee_pos_W =
                    final_com_.lin.p() + w_R_b * model_.kinematic_model_->GetNominalStanceInBase().at(ee);
            double x = final_ee_pos_W.x();
            double y = final_ee_pos_W.y();
            double z = terrain_->GetHeight(x, y);
            nodes->SetByLinearInterpolation(initial_ee_W_.at(ee), Vector3d(x, y, z), T);

            nodes->AddStartBound(kPos, {X, Y, Z}, initial_ee_W_.at(ee));
            vars.push_back(nodes);
        }

        return vars;
    }

    std::vector<NodesVariablesPhaseBased::Ptr>
    NlpFormulationMomentum::MakeForceVariables() const {
        std::vector<NodesVariablesPhaseBased::Ptr> vars;

        double T = params_.GetTotalTime();
        for (int ee = 0; ee < params_.GetEECount(); ee++) {
            auto nodes = std::make_shared<NodesVariablesEEForce>(
                    params_.GetPhaseCount(ee),
                    params_.ee_in_contact_at_start_.at(ee),
                    id::EEForceNodes(ee),
                    params_.force_polynomials_per_stance_phase_);

            // initialize with mass of robot distributed equally on all legs
            double m = model_.srbd_momentum_model_->m();
            double g = model_.srbd_momentum_model_->g();

            Vector3d f_stance(0.0, 0.0, m * g / params_.GetEECount());
            nodes->SetByLinearInterpolation(f_stance, f_stance, T); // stay constant
            vars.push_back(nodes);
        }

        return vars;
    }

    std::vector<PhaseDurations::Ptr>
    NlpFormulationMomentum::MakeContactScheduleVariables() const {
        std::vector<PhaseDurations::Ptr> vars;

        for (int ee = 0; ee < params_.GetEECount(); ee++) {
            auto var = std::make_shared<PhaseDurations>(ee,
                                                        params_.ee_phase_durations_.at(ee),
                                                        params_.ee_in_contact_at_start_.at(ee),
                                                        params_.bound_phase_duration_.first,
                                                        params_.bound_phase_duration_.second);
            vars.push_back(var);
        }

        return vars;
    }

    NlpFormulationMomentum::ConstraintPtrVec
    NlpFormulationMomentum::GetConstraints(const SplineHolderExtended &spline_holder_extended) const {
        ConstraintPtrVec constraints;
        for (auto name : params_.constraints_)
            for (auto c : GetConstraint(name, spline_holder_extended))
                constraints.push_back(c);

        return constraints;
    }

    NlpFormulationMomentum::ConstraintPtrVec
    NlpFormulationMomentum::GetConstraint(Parameters::ConstraintName name,
                                          const SplineHolderExtended &s) const {
        switch (name) {
            case Parameters::Dynamic:
                return MakeDynamicConstraintWithMomentum(s);
            case Parameters::EndeffectorRom:
                return MakeRangeOfMotionBoxConstraint(s);
            case Parameters::BaseRom:
                return MakeBaseRangeOfMotionConstraint(s);
            case Parameters::TotalTime:
                return MakeTotalTimeConstraint();
            case Parameters::Terrain:
                return MakeTerrainConstraint();
            case Parameters::Force:
                return MakeForceConstraint();
            case Parameters::Swing:
                return MakeSwingConstraint();
            case Parameters::BaseAcc:
                return MakeBaseAccConstraint(s);
            case Parameters::CoMTracking:
                return MakeCoMTrackingConstraint(s);
            default:
                throw std::runtime_error("constraint not defined!");
        }
    }

    NlpFormulationMomentum::ConstraintPtrVec
    NlpFormulationMomentum::MakeCoMTrackingConstraint(const SplineHolderExtended &s) const
    {
        return {std::make_shared<COMTrackingConstraint>(params_.target_, params_.GetTotalTime(), params_.dt_cost_com_tracking, s)};
    }


    NlpFormulationMomentum::ConstraintPtrVec
    NlpFormulationMomentum::MakeBaseRangeOfMotionConstraint(const SplineHolderExtended &s) const {
        return {std::make_shared<BaseMotionConstraintWithMomentum>(params_.GetTotalTime(),
                                                                   params_.dt_constraint_base_motion_,
                                                                   s)};
    }

    NlpFormulationMomentum::ConstraintPtrVec
    NlpFormulationMomentum::MakeDynamicConstraint(const SplineHolder &s) const {
        auto constraint = std::make_shared<DynamicConstraint>(model_.dynamic_model_,
                                                              params_.GetTotalTime(),
                                                              params_.dt_constraint_dynamic_,
                                                              s);
        return {constraint};
    }

    NlpFormulationMomentum::ConstraintPtrVec
    NlpFormulationMomentum::MakeDynamicConstraintWithMomentum(const SplineHolderExtended &s) const {
        auto constraint = std::make_shared<DynamicConstraintWithMomentum>(model_.srbd_momentum_model_,
                                                              params_.GetTotalTime(),
                                                              params_.dt_constraint_dynamic_,
                                                              s);
        return {constraint};
    }

    NlpFormulationMomentum::ConstraintPtrVec
    NlpFormulationMomentum::MakeRangeOfMotionBoxConstraint(const SplineHolderExtended &s) const {
        ConstraintPtrVec c;

        for (int ee = 0; ee < params_.GetEECount(); ee++) {
            auto rom = std::make_shared<RangeOfMotionConstraintWithMomentum>(model_.kinematic_model_,
                                                                             params_.GetTotalTime(),
                                                                             params_.dt_constraint_range_of_motion_,
                                                                             ee,
                                                                             s);
            c.push_back(rom);
        }

        return c;
    }

    NlpFormulationMomentum::ConstraintPtrVec
    NlpFormulationMomentum::MakeTotalTimeConstraint() const {
        ConstraintPtrVec c;
        double T = params_.GetTotalTime();

        for (int ee = 0; ee < params_.GetEECount(); ee++) {
            auto duration_constraint = std::make_shared<TotalDurationConstraint>(T, ee);
            c.push_back(duration_constraint);
        }

        return c;
    }

    NlpFormulationMomentum::ConstraintPtrVec
    NlpFormulationMomentum::MakeTerrainConstraint() const {
        ConstraintPtrVec constraints;

        for (int ee = 0; ee < params_.GetEECount(); ee++) {
            auto c = std::make_shared<TerrainConstraint>(terrain_, id::EEMotionNodes(ee));
            constraints.push_back(c);
        }

        return constraints;
    }

    NlpFormulationMomentum::ConstraintPtrVec
    NlpFormulationMomentum::MakeForceConstraint() const {
        ConstraintPtrVec constraints;

        for (int ee = 0; ee < params_.GetEECount(); ee++) {
            auto c = std::make_shared<ForceConstraint>(terrain_,
                                                       params_.force_limit_in_normal_direction_,
                                                       ee);
            constraints.push_back(c);
        }

        return constraints;
    }

    NlpFormulationMomentum::ConstraintPtrVec
    NlpFormulationMomentum::MakeSwingConstraint() const {
        ConstraintPtrVec constraints;

        for (int ee = 0; ee < params_.GetEECount(); ee++) {
            auto swing = std::make_shared<SwingConstraint>(id::EEMotionNodes(ee));
            constraints.push_back(swing);
        }

        return constraints;
    }

    NlpFormulationMomentum::ConstraintPtrVec
    NlpFormulationMomentum::MakeBaseAccConstraint(const SplineHolderExtended &s) const {
        ConstraintPtrVec constraints;

        constraints.push_back(std::make_shared<SplineAccConstraint>
                                      (s.com_linear_, id::base_lin_nodes));

        constraints.push_back(std::make_shared<SplineAccConstraint>
                                      (s.com_angular_, id::base_ang_nodes));

        return constraints;
    }

    NlpFormulationMomentum::ConstraintPtrVec
    NlpFormulationMomentum::GetCosts() const {
        ConstraintPtrVec costs;
        for (const auto &pair : params_.costs_)
            for (auto c : GetCost(pair.first, pair.second))
                costs.push_back(c);

        return costs;
    }

    NlpFormulationMomentum::CostPtrVec
    NlpFormulationMomentum::GetCost(const Parameters::CostName &name, double weight) const {
        switch (name) {
            case Parameters::ForcesCostID:
                return MakeForcesCost(weight);
            case Parameters::EEMotionCostID:
                return MakeEEMotionCost(weight);
            default:
                throw std::runtime_error("cost not defined!");
        }
    }

    NlpFormulationMomentum::CostPtrVec
    NlpFormulationMomentum::MakeForcesCost(double weight) const {
        CostPtrVec cost;

        for (int ee = 0; ee < params_.GetEECount(); ee++)
            cost.push_back(std::make_shared<NodeCost>(id::EEForceNodes(ee), kPos, Z, weight));

        return cost;
    }

    NlpFormulationMomentum::CostPtrVec
    NlpFormulationMomentum::MakeEEMotionCost(double weight) const {
        CostPtrVec cost;

        for (int ee = 0; ee < params_.GetEECount(); ee++) {
            cost.push_back(std::make_shared<NodeCost>(id::EEMotionNodes(ee), kVel, X, weight));
            cost.push_back(std::make_shared<NodeCost>(id::EEMotionNodes(ee), kVel, Y, weight));
        }

        return cost;
    }


    NlpFormulationMomentum::ConstraintPtrVec
    NlpFormulationMomentum::GetTrackingCosts(const SplineHolderExtended& s) const
    {
      ConstraintPtrVec costs;
      for (const auto& tuple : params_.costs_track_)
        for (auto c : GetTrackingCost(std::get<0>(tuple), std::get<1>(tuple), std::get<2>(tuple), s))
          costs.push_back(c);
    
      return costs;
    }
    
    NlpFormulationMomentum::CostPtrVec
    NlpFormulationMomentum::GetTrackingCost(const Parameters::CostName& name, double weight, std::vector<Eigen::Vector3d> target,
                                    const SplineHolderExtended& s) const
    {
      switch (name) {
    //    case Parameters::EETrackingCostID: return MakeEETrackingCost(target, weight, s);
        case Parameters::CoMTrackingCostID: return MakeCOMTrackingCost(target, weight, s);
    //    case Parameters::AMTrackingCostID: return MakeAMTrackingCost(target, weight, s);
    //    case Parameters::LMTrackingCostID: return MakeLMTrackingCost(target, weight, s);
    //    case Parameters::MoITrackingCostID: return MakeMoITrackingCost(target, weight, s);
        default: throw std::runtime_error("cost not defined!");
      }
    }


    NlpFormulationMomentum::CostPtrVec
    NlpFormulationMomentum::MakeCOMTrackingCost(std::vector<Eigen::Vector3d> target, double weight, const SplineHolderExtended& s) const{
      ConstraintPtrVec cost;

      cost.push_back(std::make_shared<COMTrackingCost>(target, weight, params_.GetTotalTime(),
                                                       params_.dt_cost_com_tracking, s));

      return cost;
    }




    NlpFormulationMomentum::ConstraintPtrVec
    NlpFormulationMomentum::GetNodeTrackingCosts() const {
      ConstraintPtrVec costs;
      for (const auto &tuple : params_.costs_node_track_)
        for (auto c : GetNodeTrackingCost(std::get<0>(tuple), std::get<1>(tuple), std::get<2>(tuple)))
          costs.push_back(c);
      return costs;
    }

    NlpFormulationMomentum::CostPtrVec
    NlpFormulationMomentum::GetNodeTrackingCost(const Parameters::CostName &name, double weight, const std::vector<Eigen::Vector3d>& target) const {
      switch (name) {
        case Parameters::CoMNodeTrackingCostID:
          return MakeCOMNodeTrackingCost(target, weight);
        default:
          throw std::runtime_error("cost not defined!");
      }
    }

    NlpFormulationMomentum::CostPtrVec
    NlpFormulationMomentum::MakeCOMNodeTrackingCost(const std::vector<Eigen::Vector3d>& target, double weight) const {
      CostPtrVec cost;

//      cost.push_back(std::make_shared<NodeTrackingCost>(id::com_lin_nodes, kPos, X, target, weight));
//      cost.push_back(std::make_shared<NodeTrackingCost>(id::com_lin_nodes, kPos, Y, target, weight));
      cost.push_back(std::make_shared<NodeTrackingCost>(id::com_lin_nodes, kPos, Z, target, weight));

      return cost;
    }


    // TODO(Bruce): Implement this member function
    /**
     * @brief take a splineHolderExtended pointer and returns a pointer to a SplineHolder pointer,
     * removes the ellipsoid variables, and angular momentum variables.
     * @param s
     * @return
     */
    SplineHolder&
    NlpFormulationMomentum::SplineHolderExtendedToRegular(const towr::SplineHolderExtended &s) const {
//        SplineHolder s_regular = SplineHolder();
//
//        return s_regular;
    }

} /* namespace towr */
