//
// Created by bruce on 3/22/21.
//

#include <cmath>
#include <iostream>

#include <towr/terrain/examples/height_map_examples.h>
#include <towr/nlp_formulation_with_momentum.h>
#include <ifopt/ipopt_solver.h>


using namespace towr;

// A minimal example how to build a trajectory optimization problem using TOWR.
//
// The more advanced example that includes ROS integration, GUI, rviz
// visualization and plotting can be found here:
// towr_ros/src/towr_ros_app.cc
int main()
{
  NlpFormulationMomentum formulation;

  // terrain
  formulation.terrain_ = std::make_shared<FlatGround>(0.0);

  // Kinematic limits and dynamic parameters of the hopper
  formulation.model_ = RobotModel(RobotModel::MonopedWithMomentum);


//  Node& debug_node = formulation.initial_am_.values;
//  debug_node = formulation.initial_am_.values;
//  std::cout << "initial am node values == " << debug_node.at(kPos) << std::endl;

  // set the initial position of the hopper
  formulation.initial_com_.lin.at(kPos).z() = 0.5;
  formulation.initial_ee_W_.push_back(Eigen::Vector3d::Zero());

  // set the initial axes and orientation of the ellipsoid
  formulation.initial_ellipsoid_state_.axis.at(kPos) << 1.0, 1.0, 1.0;
//  formulation.initial_ellipsoid_state_.axis.at(kVel) << 0.0, 0.0, 0.0;
//  formulation.initial_ellipsoid_state_.orient.at(kPos) << 0.0, 0.0, 0.0;
//  formulation.initial_ellipsoid_state_.orient.at(kVel) << 0.0, 0.0, 0.0;

  // set the initial angular momentum of the ellipsoid
//  std::cout << "initial_am_.value.at(kPos) == " << formulation.initial_am_.values.at(kVel) << std::endl;
  formulation.initial_am_.values.at(kPos).x() = 0.0;
  formulation.initial_am_.values.at(kVel) << 0.0, 0.0, 0.0;


  // define the desired goal state of the hopper
  formulation.final_com_.lin.at(towr::kPos) << 1.0, 0.0, 0.5;

  // set the final axes and orientation of the ellipsoid
  formulation.final_ellipsoid_state_.axis.at(kPos) << 1.0, 1.0, 1.0;
  formulation.final_ellipsoid_state_.axis.at(kVel) << 0.0, 0.0, 0.0;
  formulation.final_ellipsoid_state_.orient.at(kPos) << 0.0, 0.0, 0.0;
  formulation.final_ellipsoid_state_.orient.at(kVel) << 0.0, 0.0, 0.0;

  // set the final angular momentum of the ellipsoid
  formulation.final_am_.values.at(kPos) << 0.0, 0.0, 0.0;
  formulation.final_am_.values.at(kVel) << 0.0, 0.0, 0.0;

  // Parameters that define the motion. See c'tor for default values or
  // other values that can be modified.
  // First we define the initial phase durations, that can however be changed
  // by the optimizer. The number of swing and stance phases however is fixed.
  // alternating stance and swing:     ____-----_____-----_____-----_____
  formulation.params_.ee_phase_durations_.push_back({0.4, 0.2, 0.4, 0.2, 0.4, 0.2, 0.2});
  formulation.params_.ee_in_contact_at_start_.push_back(true);

  // TODO: finish initializing _am_, and _ellipsoid_state_, to avoid seg fault (done)

  // Initialize the nonlinear-programming problem with the variables,
  // constraints and costs.
  ifopt::Problem nlp;
  SplineHolderExtended solution;
  for (auto c : formulation.GetVariableSets(solution))
    nlp.AddVariableSet(c);
  for (auto c : formulation.GetConstraints(solution))
    nlp.AddConstraintSet(c);
  for (auto c : formulation.GetCosts())
    nlp.AddCostSet(c);

  // You can add your own elements to the nlp as well, simply by calling:
  // nlp.AddVariablesSet(your_custom_variables);
  // nlp.AddConstraintSet(your_custom_constraints);

  // Choose ifopt solver (IPOPT or SNOPT), set some parameters and solve.
  // solver->SetOption("derivative_test", "first-order");
  auto solver = std::make_shared<ifopt::IpoptSolver>();
  solver->SetOption("jacobian_approximation", "exact"); // "finite difference-values"
  solver->SetOption("max_cpu_time", 20.0);
  solver->SetOption("print_level", 5);
  solver->Solve(nlp);

  // Can directly view the optimization variables through:
  // Eigen::VectorXd x = nlp.GetVariableValues()
  // However, it's more convenient to access the splines constructed from these
  // variables and query their values at specific times:
  using namespace std;
  cout.precision(2);
  nlp.PrintCurrent(); // view variable-set, constraint violations, indices,...
  cout << fixed;
  cout << "\n====================\nMonoped trajectory:\n====================\n";

  double t = 0.0;
  while (t<=solution.com_linear_->GetTotalTime() + 1e-5) {
    cout << "t=" << t << "\n";
    cout << "Base linear position x,y,z:   \t";
    cout << solution.com_linear_->GetPoint(t).p().transpose() << "\t[m]" << endl;

    cout << "Base Euler roll, pitch, yaw:  \t";
    Eigen::Vector3d rad = solution.com_angular_->GetPoint(t).p();
    cout << (rad/M_PI*180).transpose() << "\t[deg]" << endl;

    cout << "Angular momentum hx,hy,hz:   \t";
    cout << solution.angular_momentum_->GetPoint(t).p().transpose() << "\t[unit]" << endl;

    cout << "Ellipsoid semi-axes ex,ey,ez:   \t";
    cout << solution.ellipsoid_axes_->GetPoint(t).p().transpose() << "\t[m]" << endl;

    cout << "Ellipsoid Euler roll, pitch, yaw:  \t";
    Eigen::Vector3d ellip_rad = solution.ellipsoid_orient_->GetPoint(t).p();
    cout << (ellip_rad/M_PI*180).transpose() << "\t[deg]" << endl;

    cout << "Foot position x,y,z:          \t";
    cout << solution.ee_motion_.at(0)->GetPoint(t).p().transpose() << "\t[m]" << endl;

    cout << "Contact force x,y,z:          \t";
    cout << solution.ee_force_.at(0)->GetPoint(t).p().transpose() << "\t[N]" << endl;

    bool contact = solution.phase_durations_.at(0)->IsContactPhase(t);
    std::string foot_in_contact = contact? "yes" : "no";
    cout << "Foot in contact:              \t" + foot_in_contact << endl;

    cout << endl;

    t += 0.2;
  }
}
