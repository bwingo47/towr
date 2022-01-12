//
// Created by bruce on 3/18/21.
//

#ifndef TOWR_SPLINE_HOLDER_EXTENDED_H
#define TOWR_SPLINE_HOLDER_EXTENDED_H

#include "phase_durations.h"
#include "node_spline.h"
#include "nodes_variables.h"
#include "nodes_variables_phase_based.h"

namespace towr {

/**
 * @brief Builds splines from node values (pos/vel) and durations.
 *
 * These splines are linked to the optimization variables, so change as the
 * nodes or durations change. This is a convenience class to not have
 * to construct the splines from the variables new every time.
 */
    struct SplineHolderExtended {
        /**
         * @brief Fully construct all splines.
         * @param base_lin  The nodes describing the base linear motion.
         * @param base_ang  The nodes describing the base angular motion.
         * @param base_poly_durations The durations of each base polynomial.
         * @param ee_motion The nodes describing the endeffector motions.
         * @param ee_force  The nodes describing the endeffector forces.
         * @param phase_durations  The phase durations of each endeffector.
         * @param ee_durations_change  True if the ee durations are optimized over.
         */

        using Ptr = std::shared_ptr<SplineHolderExtended>;

        SplineHolderExtended (NodesVariables::Ptr com_lin,
                      NodesVariables::Ptr com_ang,
                      const std::vector<double>& com_poly_durations,
                // add ellipsoid variables
                      NodesVariables::Ptr ellipsoid_axes,
                      NodesVariables::Ptr ellipsoid_orient,
                      const std::vector<double>& ellipsoid_poly_durations,
                // add angular momentum
                      NodesVariables::Ptr angular_momentum,
                      const std::vector<double>& am_poly_durations,
                //std::vector because in general more than 1 eef.
                //Each NodesVariablesPhaseBased::Ptr corresponds to
                //the ee_motion/force of a single eef.
                      std::vector<NodesVariablesPhaseBased::Ptr> ee_motion,
                      std::vector<NodesVariablesPhaseBased::Ptr> ee_force,
                      std::vector<PhaseDurations::Ptr> phase_durations,
                      bool ee_durations_change);

        SplineHolderExtended (NodesVariables::Ptr com_lin,
                              NodesVariables::Ptr com_ang,
                              const std::vector<double>& com_poly_durations,
            // add angular momentum
                              NodesVariables::Ptr angular_momentum,
                              const std::vector<double>& am_poly_durations,
            //std::vector because in general more than 1 eef.
            //Each NodesVariablesPhaseBased::Ptr corresponds to
            //the ee_motion/force of a single eef.
                              std::vector<NodesVariablesPhaseBased::Ptr> ee_motion,
                              std::vector<NodesVariablesPhaseBased::Ptr> ee_force,
                              std::vector<PhaseDurations::Ptr> phase_durations,
                              bool ee_durations_change);

        SplineHolderExtended (NodesVariables::Ptr com_lin,
                              NodesVariables::Ptr com_ang,
                              const std::vector<double>& com_poly_durations,
            // add ellipsoid variables
                              NodesVariables::Ptr ellipsoid_axes,
                              NodesVariables::Ptr ellipsoid_orient,
                              const std::vector<double>& ellipsoid_poly_durations,
            // add angular momentum phase based
                              NodesVariablesPhaseBased::Ptr angular_momentum_phase_based,
            //std::vector because in general more than 1 eef.
            //Each NodesVariablesPhaseBased::Ptr corresponds to
            //the ee_motion/force of a single eef.
                              std::vector<NodesVariablesPhaseBased::Ptr> ee_motion,
                              std::vector<NodesVariablesPhaseBased::Ptr> ee_force,
                              std::vector<PhaseDurations::Ptr> phase_durations,
                              bool ee_durations_change);

        /**
         * @brief Attention, nothing initialized.
         */
        SplineHolderExtended () = default;

        NodeSpline::Ptr com_linear_;
        NodeSpline::Ptr com_angular_;

        NodeSpline::Ptr ellipsoid_axes_;
        NodeSpline::Ptr ellipsoid_orient_;

        NodeSpline::Ptr angular_momentum_;

//        NodeSpline::Ptr angular_momentum_phase_based_;

        std::vector<NodeSpline::Ptr> ee_motion_;
        std::vector<NodeSpline::Ptr> ee_force_;
        std::vector<PhaseDurations::Ptr> phase_durations_;
    };

} /* namespace towr */


#endif //TOWR_SPLINE_HOLDER_EXTENDED_H
