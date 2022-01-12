/******************************************************************************
Copyright (c) 2018, Alexander W. Winkler. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
******************************************************************************/

#include <towr/variables/node_spline.h>

#include <towr/variables/nodes_variables.h>

namespace towr {

NodeSpline::NodeSpline(NodeSubjectPtr const node_variables,
                       const VecTimes& polynomial_durations)
    :   Spline(polynomial_durations, node_variables->GetDim()),
        NodesObserver(node_variables)
{
  UpdateNodes();
  jac_wrt_nodes_structure_ = Jacobian(node_variables->GetDim(), node_variables->GetRows());
}

void
NodeSpline::UpdateNodes ()
{
  for (int i=0; i<cubic_polys_.size(); ++i) {
    auto nodes = node_values_->GetBoundaryNodes(i);
    cubic_polys_.at(i).SetNodes(nodes.front(), nodes.back());
  }

  UpdatePolynomialCoeff();
}

int
NodeSpline::GetNodeVariablesCount() const
{
  return node_values_->GetRows();
}

NodeSpline::Jacobian
NodeSpline::GetJacobianWrtNodes (double t_global, Dx dxdt) const
{
  int id; double t_local;
  std::tie(id, t_local) = GetLocalTime(t_global, GetPolyDurations());

  return GetJacobianWrtNodes(id, t_local, dxdt);
}

NodeSpline::Jacobian
NodeSpline::GetJacobianWrtNodes (int id, double t_local, Dx dxdt) const
{
  /**
   * jac is a copy of jac_wrt_nodes_structure_, which is an empty sparse matrix.
   * This ensures that this const member function does not modify the member variable
   * jac_wrt_nodes_structure_
   */
  Jacobian jac = jac_wrt_nodes_structure_;
  FillJacobianWrtNodes(id, t_local, dxdt, jac, false);

  // needed to avoid Eigen::assert failure "wrong storage order" triggered
  // in dynamic_constraint.cc

  /**
   * compress the sparse matrix, such that no empty space left in storage.
   */
  jac.makeCompressed();

  return jac;
}

void
NodeSpline::FillJacobianWrtNodes (int poly_id, double t_local, Dx dxdt,
                                  Jacobian& jac, bool fill_with_zeros) const
{
  /**
   * DSP: derivative wrt start node's position
   * DSV: derivative wrt start node's velocity
   * DEP: derivative wrt end node's position
   * DEV: derivative wrt end node's velocity
   *
   * example shape of returned 'jac' when 'node_variables->GetDim()' == 3.
   * when idx is only associated to one nvi, i.e. for NodeVariablesAll,
   *             DSP  0   0   DSV  0   0    DEP  0   0   DEV  0   0
   * jac = [ ...  0  DSP  0    0  DSV  0     0  DEP  0    0  DEV  0  ... ]
   *              0   0  DSP   0  0   DSV    0   0  DEP   0  0   DEV
   *
   * when idx is associated to multiple nvi(s), i.e. for NodeVariablesPhaseBased,
   * say for ee motion, one optimization variable corresponds to both the start
   * and end node x-position for a constant spline. For this case,
   *             DSP + DEP     0         0     DSV + DEV     0         0
   * jac = [ ...     0     DSP + DEP     0         0     DSV + DEV     0     ... ]
   *                 0         0     DSP + DEP     0         0     DSV + DEV
   *
   */
  for (int idx=0; idx<jac.cols(); ++idx) {
    for (auto nvi : node_values_->GetNodeValuesInfo(idx)) {
      for (auto side : {NodesVariables::Side::Start, NodesVariables::Side::End}) { // every jacobian is affected by two nodes
        int node = node_values_->GetNodeId(poly_id, side);

        if (node == nvi.id_) {
          double val = 0.0;

          if (side == NodesVariables::Side::Start)
            /// DSP when nvi.deriv_ == KPoS or DSV when nvi.deriv_ == KVel
            val = cubic_polys_.at(poly_id).GetDerivativeWrtStartNode(dxdt, nvi.deriv_, t_local);
          else if (side == NodesVariables::Side::End)
            /// DEP when nvi.deriv_ == KPoS or DEV when nvi.deriv_ == KVel
            val = cubic_polys_.at(poly_id).GetDerivativeWrtEndNode(dxdt, nvi.deriv_, t_local);
          else
            assert(false); // this shouldn't happen

          // if only want structure
          if (fill_with_zeros)
            val = 0.0;

          jac.coeffRef(nvi.dim_, idx) += val;
        }
      }
    }
  }
}

} /* namespace towr */
