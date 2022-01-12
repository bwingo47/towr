//
// Created by bruce on 5/28/21.
//

#include "towr/costs/time_discretization_cost.h"

#include <cmath>
#include <iomanip>

namespace towr {


TimeDiscretizationCost::TimeDiscretizationCost (double T, double dt,
                                                std::string name)
    :CostTerm(name)
{
  double t = 0.0;
  dts_ = {t};

  for (int i=0; i<floor(T/dt); ++i) {
    t += dt;
    dts_.push_back(t);
  }

//        dts_.push_back(T); // also ensure constraints at very last node/time.
}

TimeDiscretizationCost::TimeDiscretizationCost (const VecTimes& times,
                                                std::string name)
    :CostTerm(name) // just placeholder
{
  dts_ = times;
}

int
TimeDiscretizationCost::GetNumberOfNodes () const
{
  return dts_.size();
}

//    TimeDiscretizationCost::VectorXd
//    TimeDiscretizationCost::GetValues () const
//    {
//        VectorXd cost(1);
//        cost(0) = GetCost();
//        return cost;
//    }
//
//    TimeDiscretizationCost::VecBound
//    TimeDiscretizationCost::GetBounds () const
//    {
//        return VecBound(GetRows(), ifopt::NoBound);
//    }

//    double
//    TimeDiscretizationCost::GetCost() const{
//        assert(1 == 2);
//        double cost=0;
//        return cost;
//    }

void
TimeDiscretizationCost:: FillJacobianBlock (std::string var_set,
                                            Jacobian& jac) const
{
  int k = 0;
  for (double t : dts_)
    UpdateJacobianAtInstance(t, k++, var_set, jac);
}

//    void
//    TimeDiscretizationCost::Print(double tol, int& index) const {
//        // only one scalar cost value
//        double cost = GetValues()(0);
//
//        std::cout.precision(2);
//        std::cout << std::fixed
//                  << std::left
//                  << std::setw(30) << GetName()
//                  << std::right
//                  << std::setw(4) << GetRows()
//                  << std::setw(9) << index
//                  << std::setfill ('.')
//                  << std::setw(7) << index+GetRows()-1
//                  << std::setfill (' ')
//                  << std::setw(12) << cost
//                  << std::endl;
//    }

} /* namespace towr */

