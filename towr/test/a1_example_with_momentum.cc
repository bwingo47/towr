//
// Created by bruce on 3/25/21.
//

#include <cmath>
#include <iostream>

#include <towr/terrain/examples/height_map_examples.h>
#include <towr/nlp_formulation_with_momentum.h>
#include <ifopt/ipopt_solver.h>
#include <towr/initialization/gait_generator.h>

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

  // Kinematic limits and dynamic parameters of A1
//  formulation.model_ = RobotModel(RobotModel::A1);
  formulation.model_ = RobotModel(RobotModel::A1WithMomentum);

  // set the initial position of the A1
  auto nominal_stance_B = formulation.model_.kinematic_model_->GetNominalStanceInBase();
  double z_ground = 0.0;
  formulation.initial_ee_W_ =  nominal_stance_B;
  std::for_each(formulation.initial_ee_W_.begin(), formulation.initial_ee_W_.end(),
                [&](Eigen::Vector3d& p){ p.z() = z_ground; } // feet at 0 height
  );
  formulation.initial_com_.lin.at(kPos).z() = - nominal_stance_B.front().z() + z_ground;

  // set the initial axes and orientation of the ellipsoid
//  formulation.initial_ellipsoid_state_.axis.at(kPos) << 1.0, 1.0, 1.0;

  // set the initial angular momentum of the ellipsoid
//  formulation.initial_am_.values.at(kPos).x() = 0.0;
//  formulation.initial_am_.values.at(kVel) << 0.0, 0.0, 0.0;

  // define the desired goal state of the hopper
  formulation.final_com_.lin.at(towr::kPos) << 0.0, 0.0, 0.309489;

  // set the final axes and orientation of the ellipsoid
//  formulation.final_ellipsoid_state_.axis.at(kPos) << 1.0, 1.0, 1.0;
//  formulation.final_ellipsoid_state_.axis.at(kVel) << 0.0, 0.0, 0.0;
//  formulation.final_ellipsoid_state_.orient.at(kPos) << 0.0, 0.0, 0.0;
//  formulation.final_ellipsoid_state_.orient.at(kVel) << 0.0, 0.0, 0.0;

  // set the final angular momentum of the ellipsoid
//  formulation.final_am_.values.at(kPos) << 0.0, 0.0, 0.0;
//  formulation.final_am_.values.at(kVel) << 0.0, 0.0, 0.0;

  // Parameters that define the motion. See c'tor for default values or
  // other values that can be modified.
  // First we define the initial phase durations, that can however be changed
  // by the optimizer. The number of swing and stance phases however is fixed.
  // alternating stance and swing:     ____-----_____-----_____-----_____
  Eigen::VectorXd rho(5);
  rho(0) = 1e6;

  double total_duration = 0.8;
  double time_step = 0.02;
  int knot_points;
  knot_points = floor(total_duration / time_step)+1;

  int idx = 0;

  Eigen::Vector3d com_pos_init = formulation.initial_com_.lin.at(kPos);
  Eigen::Vector3d com_pos_final = formulation.final_com_.lin.at(kPos);
  Eigen::Vector3d com_pos_mid(0.0, 0.0, 0.4);

  double com_pos_init_z = com_pos_init(2);

  double launch_phase_duration = 0.15;
  double aerial_phase_duration = 0.5;
  double landing_phase_duration = 0.15;
  double v_launch = aerial_phase_duration*9.81/2.0;

  Eigen::Vector3d com_pos_launch_mid(0.0, 0.0, com_pos_init_z - v_launch*(launch_phase_duration/2));
  Eigen::Vector3d com_pos_landing_mid(0.0, 0.0, com_pos_init_z - v_launch*(landing_phase_duration/2));

  std::vector<Eigen::Vector3d> com_target(knot_points);

//  for(auto com_t : com_target){
//
//
//
//    if(idx < floor(knot_points*0.5)){
//      double t_lerp = (idx*time_step)/(total_duration*0.5);
//      com_target[idx] = com_pos_init + t_lerp*(com_pos_mid - com_pos_init);
//    } else {
//      double t_lerp = (idx*time_step-total_duration/2)/(total_duration/2);
//      com_target[idx] = com_pos_mid + t_lerp*(com_pos_final - com_pos_mid);
//    }
//
////    com_target[idx] = com_pos_init + t_lerp*(com_pos_final - com_pos_init);
//
//    idx++;
//  }


  for(auto com_t : com_target){
    if(idx*time_step < launch_phase_duration){
      //launch phase
//      if(idx*time_step < launch_phase_duration/2){
//        double t_lerp = (idx*time_step)/(launch_phase_duration/2);
//        com_target[idx] = com_pos_init + t_lerp*(com_pos_launch_mid - com_pos_init);
//      } else {
//        double t_lerp = (idx*time_step-launch_phase_duration/2)/(launch_phase_duration/2);
//        com_target[idx] = com_pos_launch_mid + t_lerp*(com_pos_init - com_pos_launch_mid);
//      }

      double t_lerp = (idx*time_step)/launch_phase_duration;
      com_target[idx] = com_pos_init + t_lerp*(com_pos_init - com_pos_init);

    } else if(idx*time_step > (launch_phase_duration + aerial_phase_duration)){
      //landing phase
//      if((idx*time_step-launch_phase_duration-aerial_phase_duration) < landing_phase_duration/2){
//        double t_lerp = (idx*time_step - aerial_phase_duration - launch_phase_duration)/(landing_phase_duration/2);
//        com_target[idx] = com_pos_init + t_lerp*(com_pos_landing_mid - com_pos_init);
//      } else {
//        double t_lerp = (idx*time_step - aerial_phase_duration - launch_phase_duration - (landing_phase_duration/2))/(landing_phase_duration/2);
//        com_target[idx] = com_pos_landing_mid + t_lerp*(com_pos_final - com_pos_landing_mid);
//      }

      double t_lerp = (idx*time_step - aerial_phase_duration - launch_phase_duration)/landing_phase_duration;
      com_target[idx] = com_pos_init + t_lerp*(com_pos_final - com_pos_init);

    } else {
      //aerial phase, compute quadratic
      double t_aerial = (idx*time_step - launch_phase_duration);
      com_target[idx](0) = com_pos_init(0);
      com_target[idx](1) = com_pos_init(1);
      com_target[idx](2) = (-1*9.81/2)*std::pow(t_aerial, 2) + v_launch*t_aerial + com_pos_init_z;
    }
    //          cout << "com_target_z: " << x_wbd.at(i).com_pos[2] << endl;
    idx++;
  }

  // this is to make CoM tracking as constraint work
  com_target.push_back(com_target.back());


  Parameters params(com_target, rho, time_step);

  int n_ee = 4;

  // Instead of manually defining the initial durations for each foot and
  // step, for convenience we use a GaitGenerator with some predefined gaits
  // for a variety of robots (walk, trot, pace, ...).
  auto gait_gen_ = GaitGenerator::MakeGaitGenerator(n_ee);
  auto id_gait   = static_cast<GaitGenerator::Combos>(0);
  gait_gen_->SetCombo(id_gait);

  for (int ee=0; ee<n_ee; ++ee) {
//    params.ee_phase_durations_.push_back(gait_gen_->GetPhaseDurations(total_duration, ee));
    params.ee_phase_durations_.push_back({0.15, 0.5, 0.15});

//    params.ee_in_contact_at_start_.push_back(gait_gen_->IsInContactAtStart(ee));
    params.ee_in_contact_at_start_.push_back(1);
  }
  formulation.params_ = params;

//    formulation.params_.ee_phase_durations_.push_back({0.4, 0.2, 0.4, 0.2, 0.4, 0.2, 0.2});
//    formulation.params_.ee_in_contact_at_start_.push_back(true);

  // Initialize the nonlinear-programming problem with the variables,
  // constraints and costs.
  ifopt::Problem nlp;
  SplineHolderExtended solution;
  for (auto c : formulation.GetVariableSets(solution))
    nlp.AddVariableSet(c);
  for (auto c : formulation.GetConstraints(solution))
    nlp.AddConstraintSet(c);
//  for (auto c : formulation.GetCosts())
  for (auto c : formulation.GetTrackingCosts(solution))
//  for (auto c : formulation.GetNodeTrackingCosts())
    nlp.AddCostSet(c);
//    nlp.PrintCurrent();

  // You can add your own elements to the nlp as well, simply by calling:
  // nlp.AddVariablesSet(your_custom_variables);
  // nlp.AddConstraintSet(your_custom_constraints);

  // Choose ifopt solver (IPOPT or SNOPT), set some parameters and solve.

  auto solver = std::make_shared<ifopt::IpoptSolver>();
  solver->SetOption("jacobian_approximation", "exact"); // "finite-difference-values" / "exact"
//  solver->SetOption("derivative_test", "second-order");
  solver->SetOption("max_cpu_time", 7200.0);
  solver->SetOption("print_level", 5);
  solver->SetOption("max_iter", 500);
  solver->Solve(nlp);

  // Can directly view the optimization variables through:
  // Eigen::VectorXd x = nlp.GetVariableValues()
  // However, it's more convenient to access the splines constructed from these
  // variables and query their values at specific times:
  using namespace std;
  cout.precision(2);
  nlp.PrintCurrent(); // view variable-set, constraint violations, indices,...
  cout << fixed;
  cout << "\n====================\nA1 trajectory:\n====================\n";

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

    t += 0.1;
  }
}


