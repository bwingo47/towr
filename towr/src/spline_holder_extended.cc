//
// Created by bruce on 3/18/21.
//

#include <towr/variables/spline_holder_extended.h>
#include <towr/variables/phase_spline.h>

namespace towr{

    SplineHolderExtended::SplineHolderExtended (NodesVariables::Ptr com_lin_nodes,
                                NodesVariables::Ptr com_ang_nodes,
                                const std::vector<double>& com_poly_durations,
                                NodesVariables::Ptr ellipsoid_axes_nodes,
                                NodesVariables::Ptr ellipsoid_orient_nodes,
                                const std::vector<double>& ellipsoid_poly_durations,
                                NodesVariables::Ptr am_nodes,
                                const std::vector<double>& am_poly_durations,
                                std::vector<NodesVariablesPhaseBased::Ptr> ee_motion_nodes,
                                std::vector<NodesVariablesPhaseBased::Ptr> ee_force_nodes,
                                std::vector<PhaseDurations::Ptr> phase_durations,
                                bool durations_change)
    {
        com_linear_  = std::make_shared<NodeSpline>(com_lin_nodes.get(), com_poly_durations);
        com_angular_ = std::make_shared<NodeSpline>(com_ang_nodes.get(), com_poly_durations);

        ellipsoid_axes_ = std::make_shared<NodeSpline>(ellipsoid_axes_nodes.get(), ellipsoid_poly_durations);
        ellipsoid_orient_ = std::make_shared<NodeSpline>(ellipsoid_orient_nodes.get(), ellipsoid_poly_durations);

        angular_momentum_ = std::make_shared<NodeSpline>(am_nodes.get(), am_poly_durations);

        phase_durations_ = phase_durations;

        for (uint ee=0; ee<ee_motion_nodes.size(); ++ee) {
            if (durations_change) {
                // spline that changes the polynomial durations (affects Jacobian)
                ee_motion_.push_back(std::make_shared<PhaseSpline>(ee_motion_nodes.at(ee), phase_durations.at(ee).get()));
                ee_force_.push_back(std::make_shared<PhaseSpline>(ee_force_nodes.at(ee), phase_durations.at(ee).get()));
            } else {
                // spline without changing the polynomial durations
                auto ee_motion_poly_durations = ee_motion_nodes.at(ee)->ConvertPhaseToPolyDurations(phase_durations.at(ee)->GetPhaseDurations());
                auto ee_force_poly_durations = ee_force_nodes.at(ee)->ConvertPhaseToPolyDurations(phase_durations.at(ee)->GetPhaseDurations());

                ee_motion_.push_back(std::make_shared<NodeSpline>(ee_motion_nodes.at(ee).get(), ee_motion_poly_durations));
                ee_force_.push_back (std::make_shared<NodeSpline>(ee_force_nodes.at(ee).get(), ee_force_poly_durations));
            }
        }
    }



    SplineHolderExtended::SplineHolderExtended (NodesVariables::Ptr com_lin_nodes,
                                                NodesVariables::Ptr com_ang_nodes,
                                                const std::vector<double>& com_poly_durations,
                                                NodesVariables::Ptr am_nodes,
                                                const std::vector<double>& am_poly_durations,
                                                std::vector<NodesVariablesPhaseBased::Ptr> ee_motion_nodes,
                                                std::vector<NodesVariablesPhaseBased::Ptr> ee_force_nodes,
                                                std::vector<PhaseDurations::Ptr> phase_durations,
                                                bool durations_change)
    {
        com_linear_  = std::make_shared<NodeSpline>(com_lin_nodes.get(), com_poly_durations);
        com_angular_ = std::make_shared<NodeSpline>(com_ang_nodes.get(), com_poly_durations);

        angular_momentum_ = std::make_shared<NodeSpline>(am_nodes.get(), am_poly_durations);

        phase_durations_ = phase_durations;

        for (uint ee=0; ee<ee_motion_nodes.size(); ++ee) {
            if (durations_change) {
                // spline that changes the polynomial durations (affects Jacobian)
                ee_motion_.push_back(std::make_shared<PhaseSpline>(ee_motion_nodes.at(ee), phase_durations.at(ee).get()));
                ee_force_.push_back(std::make_shared<PhaseSpline>(ee_force_nodes.at(ee), phase_durations.at(ee).get()));
            } else {
                // spline without changing the polynomial durations
                auto ee_motion_poly_durations = ee_motion_nodes.at(ee)->ConvertPhaseToPolyDurations(phase_durations.at(ee)->GetPhaseDurations());
                auto ee_force_poly_durations = ee_force_nodes.at(ee)->ConvertPhaseToPolyDurations(phase_durations.at(ee)->GetPhaseDurations());

                ee_motion_.push_back(std::make_shared<NodeSpline>(ee_motion_nodes.at(ee).get(), ee_motion_poly_durations));
                ee_force_.push_back (std::make_shared<NodeSpline>(ee_force_nodes.at(ee).get(), ee_force_poly_durations));
            }
        }
    }


    SplineHolderExtended::SplineHolderExtended (NodesVariables::Ptr com_lin_nodes,
                                                NodesVariables::Ptr com_ang_nodes,
                                                const std::vector<double>& com_poly_durations,
                                                NodesVariables::Ptr ellipsoid_axes_nodes,
                                                NodesVariables::Ptr ellipsoid_orient_nodes,
                                                const std::vector<double>& ellipsoid_poly_durations,
                                                NodesVariablesPhaseBased::Ptr am_nodes_phase_based,
                                                std::vector<NodesVariablesPhaseBased::Ptr> ee_motion_nodes,
                                                std::vector<NodesVariablesPhaseBased::Ptr> ee_force_nodes,
                                                std::vector<PhaseDurations::Ptr> phase_durations,
                                                bool durations_change)
    {
        com_linear_  = std::make_shared<NodeSpline>(com_lin_nodes.get(), com_poly_durations);
        com_angular_ = std::make_shared<NodeSpline>(com_ang_nodes.get(), com_poly_durations);

        ellipsoid_axes_ = std::make_shared<NodeSpline>(ellipsoid_axes_nodes.get(), ellipsoid_poly_durations);
        ellipsoid_orient_ = std::make_shared<NodeSpline>(ellipsoid_orient_nodes.get(), ellipsoid_poly_durations);

        phase_durations_ = phase_durations;

        auto am_phase_based_poly_durations = am_nodes_phase_based->ConvertPhaseToPolyDurations(phase_durations.at(0)->GetPhaseDurations());
        angular_momentum_ = std::make_shared<NodeSpline>(am_nodes_phase_based.get(), am_phase_based_poly_durations);

        for (uint ee=0; ee<ee_motion_nodes.size(); ++ee) {
            if (durations_change) {
                // spline that changes the polynomial durations (affects Jacobian)
                ee_motion_.push_back(std::make_shared<PhaseSpline>(ee_motion_nodes.at(ee), phase_durations.at(ee).get()));
                ee_force_.push_back(std::make_shared<PhaseSpline>(ee_force_nodes.at(ee), phase_durations.at(ee).get()));
            } else {
                // spline without changing the polynomial durations
                auto ee_motion_poly_durations = ee_motion_nodes.at(ee)->ConvertPhaseToPolyDurations(phase_durations.at(ee)->GetPhaseDurations());
                auto ee_force_poly_durations = ee_force_nodes.at(ee)->ConvertPhaseToPolyDurations(phase_durations.at(ee)->GetPhaseDurations());

                ee_motion_.push_back(std::make_shared<NodeSpline>(ee_motion_nodes.at(ee).get(), ee_motion_poly_durations));
                ee_force_.push_back (std::make_shared<NodeSpline>(ee_force_nodes.at(ee).get(), ee_force_poly_durations));
            }
        }
    }

} /* namespace towr */
