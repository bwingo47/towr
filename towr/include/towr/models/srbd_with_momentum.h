//
// Created by bruce on 3/18/21.
//

#ifndef TOWR_SRBD_WITH_MOMENTUM_H
#define TOWR_SRBD_WITH_MOMENTUM_H

#include <memory>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <towr/variables/euler_converter.h>

namespace towr {

/**
 * @brief Dynamics model relating forces to base accelerations.
 *
 * This class implements a Single Rigid Body dynamics model, a reduced
 * dimensional model, lying in terms of accuracy between a Linear
 * Inverted Pendulum model and a full Centroidal or Rigid-body-dynamics model.
 *
 * This model makes the assumption that the motion of the limbs does not
 * incur significant momentum and can therefore be neglected. This eliminates
 * the nonlinear dependency on joint angles and allows to express
 * all quantities in Cartesian space.
 *
 * For the derivation and all assumptions of this model, see:
 * https://doi.org/10.3929/ethz-b-000272432
 *
 * \sa https://en.wikipedia.org/wiki/Newton%E2%80%93Euler_equations
 *
 * @ingroup Robots
 */
    class SingleRigidBodyDynamicsWithMomentum {
    public:
        using Ptr      = std::shared_ptr<SingleRigidBodyDynamicsWithMomentum>;
        using Vector3d = Eigen::Vector3d;
        using Matrix3d = Eigen::Matrix3d;
        using ComPos   = Eigen::Vector3d;
        using AngVel   = Eigen::Vector3d;
        using BaseAcc  = Eigen::Matrix<double,6,1>;
        using BaseAccWithAM = Eigen::Matrix<double,9,1>;
        using Jac      = Eigen::SparseMatrix<double, Eigen::RowMajor>;
        using EEPos    = std::vector<Eigen::Vector3d>;
        using EELoad   = EEPos;
        using EE       = uint;
        /**
         * @brief Constructs a specific model.
         * @param mass         The mass of the robot.
         * @param ee_count     The number of endeffectors/forces.
         * @param inertia_b    The elements of the 3x3 Inertia matrix around the CoM.
         *                     This matrix maps angular accelerations expressed in
         *                     base frame to moments in base frame.
         */
        SingleRigidBodyDynamicsWithMomentum (double mass, const Eigen::Matrix3d& inertia_b, int ee_count);

        /**
         * @brief Constructs a specific model.
         * @param mass      Mass of the robot.
         * @param I..       Elements of the 3x3 Inertia matrix
         * @param ee_count  Number of endeffectors/forces.
         */
        SingleRigidBodyDynamicsWithMomentum (double mass,
                                 double Ixx, double Iyy, double Izz,
                                 double Ixy, double Ixz, double Iyz,
                                 int ee_count);

        virtual ~SingleRigidBodyDynamicsWithMomentum () = default;



        /**
         * @brief Sets the current state and input of the system.
         * @param com_W        Current Center-of-Mass (x,y,z) position in world frame.
         * @param com_acc_W    Current Center-of-Mass (x,y,z) acceleration in world.
         * @param w_R_b        Current rotation from base to world frame.
         * @param omega_W      Current angular velocity in world frame.
         * @param omega_dot_W  Current angular acceleration in world frame.
         * @param force_W      Force at each foot expressed in world frame.
         * @param pos_W        Position of each foot expressed in world frame
         */
        void SetCurrent(const ComPos& com_W, const Vector3d com_acc_W,
                        const Matrix3d& w_R_b, const AngVel& omega_W, const Vector3d& omega_dot_W,
                        const EELoad& force_W, const EEPos& pos_W, const Vector3d ellipsoid_axes,
                        const Matrix3d ellipsoid_w_R_b, const Vector3d am, const Vector3d am_dot);

        /**
         * @brief  The violation of the system dynamics incurred by the current values.
         * @return The 9-dimension generalized force violation (angular + linear + am).
         */
        BaseAccWithAM GetDynamicViolation() const;

        /**
         * @brief How the base position affects the dynamic violation.
         * @param jac_base_lin_pos  The 3xn Jacobian of the base linear position.
         * @param jac_base_lin_acc  The 3xn Jacobian of the base linear acceleration.
         *
         * @return The 9xn Jacobian of dynamic violations with respect to
         *         variables defining the base linear spline (e.g. node values).
         */
        Jac GetJacobianWrtBaseLin(const Jac& jac_base_lin_pos,
                                  const Jac& jac_acc_base_lin) const;

        /**
         * @brief How the base orientation affects the dynamic violation.
         * @param base_angular  provides Euler angles Jacobians.
         * @param t  Time at which euler angles values are queried.
         *
         * @return The 9xn Jacobian of dynamic violations with respect to
         *         variables defining the base angular spline (e.g. node values).
         */
        Jac GetJacobianWrtBaseAng(const EulerConverter& base_angular,
                                  double t) const;

        /**
         * @brief How the angular momentum nodes affects the dynamic violation.
         * @param jac_angular_momentum_pos  The 3xn Jacobian of the angular momentum.
         * @param jac_angular_momentum_vel  The 3xn Jacobian of the angular momentum rate.
         *
         * @return The 9xn Jacobian of dynamic violations with respect to
         *         variables defining the angular momentum spline (e.g. node values).
         */
        Jac GetJacobianWrtAM(const Jac& jac_angular_momentum_pos,
                                  const Jac& jac_angular_momentum_vel) const;


        /**
         * @brief How the ellipsoid principle semi axes affects the dynamic violation.
         * @param jac_ellipsoid_axes_pos  The 3xn Jacobian of the ellipsoid_axes.
         *
         * @return The 9xn Jacobian of dynamic violations with respect to
         *         variables defining the ellipsoid principle semi axes spline (e.g. node values).
         */
        Jac GetJacobianWrtEllipsoidAxes(const Jac& jac_ellipsoid_axes_pos) const;


        /**
         * @brief How the ellipsoid orientation affects the dynamic violation.
         * @param ellipsoid_orient  provides Euler angles Jacobians.
         * @param t  Time at which euler angles values are queried.
         *
         * @return The 9xn Jacobian of dynamic violations with respect to
         *         variables defining the ellipsoid orientation spline (e.g. node values).
         */
        Jac GetJacobianWrtEllipsoidOrient(const EulerConverter& ellipsoid_orient,
                                          double t) const;


        /**
         * @brief How the endeffector forces affect the dynamic violation.
         * @param ee_force  The 3xn Jacobian of the foot force x,y,z.
         * @param ee        The endeffector for which the senstivity is required.
         *
         * @return The 9xn Jacobian of dynamic violations with respect to
         *         variables defining the endeffector forces (e.g. node values).
         */
        Jac GetJacobianWrtForce(const Jac& jac_force, EE) const;

        /**
         * @brief How the endeffector positions affect the dynamic violation.
         * @param ee_force  The 3xn Jacobian of the foot position x,y,z.
         * @param ee        The endeffector for which the senstivity is required.
         *
         * @return The 9xn Jacobian of dynamic violations with respect to
         *         variables defining the foot positions (e.g. node values).
         */
        Jac GetJacobianWrtEEPos(const Jac& jac_ee_pos, EE) const;

        /**
         * @returns The gravity acceleration [m/s^2] (positive)
         */
        double g() const { return g_; };

        /**
         * @returns The mass of the robot [kg].
         */
        double m() const { return m_; };

        /**
         * @brief the number of endeffectors that this robot has.
         */
        int GetEECount() const { return ee_pos_.size(); };

        // dynamics helper functions
        static Matrix3d BuildInertiaTensor( double Ixx, double Iyy, double Izz,
                                            double Ixy, double Ixz, double Iyz);

        // builds a cross product matrix out of "in", so in x v = X(in)*v
        Jac Cross(const Eigen::Vector3d& in) const;

    protected:
        ComPos com_pos_;   ///< x-y-z position of the Center-of-Mass.
        Vector3d com_acc_; ///< x-y-z acceleration of the Center-of-Mass.

        Matrix3d w_R_b_;     ///< rotation matrix from base (b) to world (w) frame.
        AngVel omega_;       ///< angular velocity expressed in world frame.
        Vector3d omega_dot_; ///< angular acceleration expressed in world frame.

        Vector3d ellipsoid_axes_;   ///< ellipsoid principle semi-axes, the 'e' vector
        Matrix3d ellipsoid_w_R_b_;  ///< rotation matrix R(v) from base coord. to world coord.
        Eigen::SparseMatrix<double, Eigen::RowMajor> I_ellipsoid;       ///< moment of inertia of the ellipsoid in world frame.
        Eigen::SparseMatrix<double, Eigen::RowMajor> I_ellipsoid_body;  ///< moment of inertia of the ellipsoid in body frame.

        Vector3d am_;               /// angular momentum vector H
        Vector3d am_dot_;           /// angular momentum rate vector H_dot

        EEPos  ee_pos_;   ///< The x-y-z position of each endeffector.
        EELoad ee_force_; ///< The endeffector force expressed in world frame.

    private:
        /** Inertia of entire robot around the CoM expressed in a frame anchored
         *  in the base.
         */
        Eigen::SparseMatrix<double, Eigen::RowMajor> I_b;

        // update ellipsoid inertia
//        Matrix3d BuildEllispoidInertiaTensor(Vector3d& ellipsoid_axes, Matrix3d& ellipsoid_w_R_b);
        void UpdateEllispoidInertiaTensor();

        double g_; ///< gravity acceleration [m/s^2]
        double m_; ///< mass of the robot
    };


} /* namespace towr */

#endif //TOWR_SRBD_WITH_MOMENTUM_H
