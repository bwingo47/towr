//
// Created by bruce on 3/18/21.
//

#include <towr/models/srbd_with_momentum.h>
#include <towr/variables/cartesian_dimensions.h>

#include <cmath>

namespace towr {

SingleRigidBodyDynamicsWithMomentum::SingleRigidBodyDynamicsWithMomentum (double mass,
                                                  double Ixx, double Iyy, double Izz,
                                                  double Ixy, double Ixz, double Iyz,
                                                  int ee_count)
    : SingleRigidBodyDynamicsWithMomentum(mass,
                              BuildInertiaTensor(Ixx, Iyy, Izz, Ixy, Ixz, Iyz),
                              ee_count)
{
}

SingleRigidBodyDynamicsWithMomentum::SingleRigidBodyDynamicsWithMomentum (double mass, const Eigen::Matrix3d& inertia_b,
                                                  int ee_count)
{
  I_b = inertia_b.sparseView();

  m_ = mass;
  g_ = 9.80665;

  com_pos_.setZero();
  com_acc_.setZero();

  w_R_b_.setIdentity();
  omega_.setZero();
  omega_dot_ .setZero();

  // TODO: does it make sense to initialize principle axes as the one vector? Should we initialize axes using inertia_b?
  ellipsoid_axes_.setOnes();
  ellipsoid_w_R_b_.setIdentity();
  UpdateEllispoidInertiaTensor();

  am_.setZero();
  am_dot_.setOnes();

  ee_force_ = EELoad(ee_count);
  ee_pos_ = EEPos(ee_count);
}


// some Eigen helper functions
SingleRigidBodyDynamicsWithMomentum::Matrix3d
SingleRigidBodyDynamicsWithMomentum::BuildInertiaTensor( double Ixx, double Iyy, double Izz,
                                                         double Ixy, double Ixz, double Iyz)
{
  Eigen::Matrix3d I;
  I <<  Ixx, -Ixy, -Ixz,
       -Ixy,  Iyy, -Iyz,
       -Ixz, -Iyz,  Izz;
  return I;
}


// TODO: double check the ellipsoid inertia formula. (checked)
void
SingleRigidBodyDynamicsWithMomentum::UpdateEllispoidInertiaTensor()
{
  double e1 = ellipsoid_axes_(0);
  double e2 = ellipsoid_axes_(1);
  double e3 = ellipsoid_axes_(2);

  double I_xx = (1.0/5.0)*m_*(std::pow(e2, 2) + std::pow(e3, 2));
  double I_yy = (1.0/5.0)*m_*(std::pow(e3, 2) + std::pow(e1, 2));
  double I_zz = (1.0/5.0)*m_*(std::pow(e1, 2) + std::pow(e2, 2));

  Eigen::Matrix3d Ib, Iw;
  Ib << I_xx,  0,     0,
         0,   I_yy,   0,
         0,    0,    I_zz;

  Iw = ellipsoid_w_R_b_*Ib*ellipsoid_w_R_b_.transpose();
  I_ellipsoid_body = Ib.sparseView();
  I_ellipsoid = Iw.sparseView();
}


// builds a cross product matrix out of "in", so in x v = X(in)*v
SingleRigidBodyDynamicsWithMomentum::Jac
SingleRigidBodyDynamicsWithMomentum::Cross(const Eigen::Vector3d& in) const
{
  SingleRigidBodyDynamicsWithMomentum::Jac out(3,3);

  out.coeffRef(0,1) = -in(2); out.coeffRef(0,2) =  in(1);
  out.coeffRef(1,0) =  in(2); out.coeffRef(1,2) = -in(0);
  out.coeffRef(2,0) = -in(1); out.coeffRef(2,1) =  in(0);

  return out;
}


void
SingleRigidBodyDynamicsWithMomentum::SetCurrent (const ComPos& com_W, const Vector3d com_acc_W,
                          const Matrix3d& w_R_b, const AngVel& omega_W, const Vector3d& omega_dot_W,
                          const EELoad& force_W, const EEPos& pos_W, const Vector3d ellipsoid_axes,
                          const Matrix3d ellipsoid_w_R_b, const Vector3d am, const Vector3d am_dot)
{
  com_pos_   = com_W;
  com_acc_   = com_acc_W;

  w_R_b_     = w_R_b;
  omega_     = omega_W;
  omega_dot_ = omega_dot_W;

  ellipsoid_axes_ = ellipsoid_axes;
  ellipsoid_w_R_b_ = ellipsoid_w_R_b;
  UpdateEllispoidInertiaTensor();

  am_ = am;
  am_dot_ = am_dot;

  ee_force_  = force_W;
  ee_pos_    = pos_W;
}

// TODO: check this implementation (checked)
SingleRigidBodyDynamicsWithMomentum::BaseAccWithAM
SingleRigidBodyDynamicsWithMomentum::GetDynamicViolation () const
{
  // https://en.wikipedia.org/wiki/Newton%E2%80%93Euler_equations

  Vector3d f_sum, tau_sum;
  f_sum.setZero(); tau_sum.setZero();

  for (int ee=0; ee<ee_pos_.size(); ++ee) {
    Vector3d f = ee_force_.at(ee);
    tau_sum += f.cross(com_pos_ - ee_pos_.at(ee));
    f_sum   += f;
  }

//  // express inertia matrix in world frame based on current body orientation
//  Jac I_w = w_R_b_.sparseView() * I_b * w_R_b_.transpose().sparseView();

  BaseAccWithAM acc;
  acc.segment(aX, k3D) = am_dot_
      - tau_sum;
  acc.segment(lX, k3D) = m()*com_acc_
      - f_sum
      - Vector3d(0.0, 0.0, -m()*g()); // gravity force
  acc.segment(hX, k3D) = I_ellipsoid*omega_
      - am_;
  return acc;
}

// TODO: check this implementation (checked)
SingleRigidBodyDynamicsWithMomentum::Jac
SingleRigidBodyDynamicsWithMomentum::GetJacobianWrtBaseLin (const Jac& jac_pos_base_lin,
                                                const Jac& jac_acc_base_lin) const
{
  // build the com jacobian
  int n = jac_pos_base_lin.cols();

  Jac jac_tau_sum(k3D, n);
  for (const Vector3d& f : ee_force_) {
    Jac jac_tau = Cross(f)*jac_pos_base_lin;
    jac_tau_sum += jac_tau;
  }

  Jac jac(k9D, n);

  /// derivative of first constraint eqn wrt com_lin_ nodes
  jac.middleRows(aX, k3D) = -jac_tau_sum;

  /// derivative of second constraint eqn wrt com_lin_ nodes
  jac.middleRows(lX, k3D) = m()*jac_acc_base_lin;
//  Jac jac_zeros(k3D, n);
//  jac_zeros.setZero();
//  jac.middleRows(hX, k3D) = jac_zeros;
  /// derivative of third constraint eqn wrt com_lin_ nodes is 0 == 0

  return jac;
}

// TODO: check this implementation
SingleRigidBodyDynamicsWithMomentum::Jac
SingleRigidBodyDynamicsWithMomentum::GetJacobianWrtBaseAng (const EulerConverter& base_euler,
                                                            double t) const
{
//  Jac I_w = w_R_b_.sparseView() * I_b * w_R_b_.transpose().sparseView();
//
//  // Derivative of R*I_b*R^T * wd
//  // 1st term of product rule (derivative of R)
//  // (partial[w_R_b_]/partial[com_angular_nodes])*I_b*(w_R_b_^{T})*(omega_dot_)
//  Vector3d v11 = I_b*w_R_b_.transpose()*omega_dot_;
//  Jac jac11 = base_euler.DerivOfRotVecMult(t, v11, false);
//
//  // 2nd term of product rule (derivative of R^T)
//  // w_R_b_*I_b*((partial[w_R_b_]/partial[com_angular_nodes])^{T})*(omega_dot_)
//  Jac jac12 = w_R_b_.sparseView()*I_b*base_euler.DerivOfRotVecMult(t, omega_dot_, true);
//
//  // 3rd term of product rule (derivative of wd)
//  // I_w*(partial[omega_dot_]/partial[com_angular_nodes])
//  Jac jac_ang_acc = base_euler.GetDerivOfAngAccWrtEulerNodes(t);
//  Jac jac13 = I_w * jac_ang_acc;
//  Jac jac1 = jac11 + jac12 + jac13;
//
//
//  // Derivative of w x Iw
//  // w x d_dn(R*I_b*R^T*w) -(I*w x d_dnw)
//  // right derivative same as above, just with velocity instead acceleration
//  Vector3d v21 = I_b*w_R_b_.transpose()*omega_;
//  Jac jac21 = base_euler.DerivOfRotVecMult(t, v21, false);
//
//  // 2nd term of product rule (derivative of R^T)
//  Jac jac22 = w_R_b_.sparseView()*I_b*base_euler.DerivOfRotVecMult(t, omega_, true);
//
//  // 3rd term of product rule (derivative of omega)
//  Jac jac_ang_vel = base_euler.GetDerivOfAngVelWrtEulerNodes(t);
//  Jac jac23 = I_w * jac_ang_vel;
//
//  Jac jac2 = Cross(omega_)*(jac21+jac22+jac23) - Cross(I_w*omega_)*jac_ang_vel;
//
//
//  // Combine the two to get sensitivity to I_w*w + w x (I_w*w)
//  int n = jac_ang_vel.cols();
//  Jac jac(k6D, n);
//  jac.middleRows(AX, k3D) = jac1 + jac2;

  // partial[omega_]/partial[com_angular_nodes]
  Jac jac_ang_vel = base_euler.GetDerivOfAngVelWrtEulerNodes(t);

  int n = jac_ang_vel.cols();
  Jac jac(k9D, n);
  /// derivative of first constraint eqn wrt com_ang_ nodes is 0 == 0
  /// derivative of second constraint eqn wrt com_ang_ nodes is 0 == 0

  /// derivative of third constraint eqn wrt com_lin_ nodes
  jac.middleRows(hX, k3D) = I_ellipsoid*jac_ang_vel;

  return jac;
}

// TODO: check this implementation (checked)
SingleRigidBodyDynamicsWithMomentum::Jac
SingleRigidBodyDynamicsWithMomentum::GetJacobianWrtForce (const Jac& jac_force, EE ee) const
{
  Vector3d r = com_pos_ - ee_pos_.at(ee);
  Jac jac_tau = -Cross(r)*jac_force;

  int n = jac_force.cols();
  Jac jac(k9D, n);

  /// derivative of first constraint eqn wrt ee_forces_nodes
  jac.middleRows(aX, k3D) = -jac_tau;

  /// derivative of second constraint eqn wrt ee_forces_nodes
  jac.middleRows(lX, k3D) = -jac_force;

  /// derivative of third constraint eqn wrt ee_forces_nodes is 0 == 0

  return jac;
}

// TODO: check this implementation (checked)
SingleRigidBodyDynamicsWithMomentum::Jac
SingleRigidBodyDynamicsWithMomentum::GetJacobianWrtEEPos (const Jac& jac_ee_pos, EE ee) const
{
  Vector3d f = ee_force_.at(ee);
  Jac jac_tau = Cross(f)*(-jac_ee_pos);

  Jac jac(k9D, jac_tau.cols());
  /// derivative of first constraint eqn wrt ee_motion_nodes
  jac.middleRows(aX, k3D) = -jac_tau;

  /// derivative of second constraint eqn wrt ee_motion_nodes is 0 == 0

  /// derivative of third constraint eqn wrt ee_motion_nodes is 0 == 0

  // linear dynamics don't depend on endeffector position.
  return jac;
}

// TODO: check this implement (checked)
SingleRigidBodyDynamicsWithMomentum::Jac
SingleRigidBodyDynamicsWithMomentum::GetJacobianWrtAM(const towr::SingleRigidBodyDynamicsWithMomentum::Jac &jac_angular_momentum_pos,
                                                      const towr::SingleRigidBodyDynamicsWithMomentum::Jac &jac_angular_momentum_vel) const
{
  int n = jac_angular_momentum_pos.cols();
  Jac jac(k9D, n);

  /// derivative of first constraint eqn wrt angular_momentum_nodes
  jac.middleRows(aX, k3D) = jac_angular_momentum_vel;

  /// derivative of second constraint eqn wrt angular_momentum_nodes is 0 == 0

  /// derivative of third constraint eqn wrt angular_momentum_nodes
  jac.middleRows(hX, k3D) = -jac_angular_momentum_pos;

  return jac;
};

// TODO: check this implementation
SingleRigidBodyDynamicsWithMomentum::Jac
SingleRigidBodyDynamicsWithMomentum::GetJacobianWrtEllipsoidAxes(const towr::SingleRigidBodyDynamicsWithMomentum::Jac &jac_ellipsoid_axes_pos) const
{
  int n = jac_ellipsoid_axes_pos.cols();
  Jac jac(k9D, n);
  // p = ellipoid_w_R_b.transpose() times am_dot
  // get the jacobian from the element-wise matrix multiplications
  Vector3d  p = ellipsoid_w_R_b_.transpose()*omega_;
  Jac jac1(3, 3);
  jac1.coeffRef(0, 0) = 0.4*m()*(ellipsoid_w_R_b_.coeffRef(0,1)*p(1)+ellipsoid_w_R_b_.coeffRef(0,2)*p(2))*ellipsoid_axes_(0);
  jac1.coeffRef(0, 1) = 0.4*m()*(ellipsoid_w_R_b_.coeffRef(0,0)*p(0)+ellipsoid_w_R_b_.coeffRef(0,2)*p(2))*ellipsoid_axes_(1);
  jac1.coeffRef(0, 2) = 0.4*m()*(ellipsoid_w_R_b_.coeffRef(0,0)*p(0)+ellipsoid_w_R_b_.coeffRef(0,1)*p(1))*ellipsoid_axes_(2);
  jac1.coeffRef(1, 0) = 0.4*m()*(ellipsoid_w_R_b_.coeffRef(1,1)*p(1)+ellipsoid_w_R_b_.coeffRef(1,2)*p(2))*ellipsoid_axes_(0);
  jac1.coeffRef(1, 1) = 0.4*m()*(ellipsoid_w_R_b_.coeffRef(1,0)*p(0)+ellipsoid_w_R_b_.coeffRef(1,2)*p(2))*ellipsoid_axes_(1);
  jac1.coeffRef(1, 2) = 0.4*m()*(ellipsoid_w_R_b_.coeffRef(1,0)*p(0)+ellipsoid_w_R_b_.coeffRef(1,1)*p(1))*ellipsoid_axes_(2);
  jac1.coeffRef(2, 0) = 0.4*m()*(ellipsoid_w_R_b_.coeffRef(2,1)*p(1)+ellipsoid_w_R_b_.coeffRef(2,2)*p(2))*ellipsoid_axes_(0);
  jac1.coeffRef(2, 1) = 0.4*m()*(ellipsoid_w_R_b_.coeffRef(2,0)*p(0)+ellipsoid_w_R_b_.coeffRef(2,2)*p(2))*ellipsoid_axes_(1);
  jac1.coeffRef(2, 2) = 0.4*m()*(ellipsoid_w_R_b_.coeffRef(2,0)*p(0)+ellipsoid_w_R_b_.coeffRef(2,1)*p(1))*ellipsoid_axes_(2);
//  jac = jac1 * jac_ellipsoid_axes_pos;
  jac.middleRows(hX, k3D) = jac1 * jac_ellipsoid_axes_pos;
  return jac;
}

// TODO: check this implement (checked)
SingleRigidBodyDynamicsWithMomentum::Jac
SingleRigidBodyDynamicsWithMomentum::GetJacobianWrtEllipsoidOrient(const towr::EulerConverter &ellipsoid_orient,
                                                                   double t) const
{
  // Derivative of R*I_b*R^T * wd
  // 1st term of product rule (derivative of R)
  // (partial[ellipsoid_w_R_b_]/partial[ellipsoid_orient_nodes])*I_ellipsoid_body*(ellipsoid_w_R_b_^{T})*(omega_)
  Vector3d v1 = I_ellipsoid_body*ellipsoid_w_R_b_.transpose()*omega_;
  Jac jac1 = ellipsoid_orient.DerivOfRotVecMult(t, v1, false);

  // 2nd term of product rule (derivative of R^T)
  // ellipsoid_w_R_b_*I_ellipsoid_body*((partial[ellipsoid_w_R_b_]/partial[ellipsoid_orient_nodes])^{T})*(omega_)
  Jac jac2 = ellipsoid_w_R_b_.sparseView()*I_ellipsoid_body*ellipsoid_orient.DerivOfRotVecMult(t, omega_, true);

  Jac jac_ang_vel = ellipsoid_orient.GetDerivOfAngVelWrtEulerNodes(t);
  int n = jac_ang_vel.cols();

  Jac jac(k9D, n);

  /// derivative of first constraint eqn wrt ellipsoid_orientation_nodes is 0 == 0

  /// derivative of second constraint eqn wrt ellipsoid_orientation_nodes is 0 == 0

  /// derivative of third constraint eqn wrt ellipsoid_orientation_nodes
  jac.middleRows(hX, k3D) = jac1 + jac2;

  return jac;
}



} /* namespace towr */


