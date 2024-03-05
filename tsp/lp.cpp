#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "ortools/linear_solver/linear_solver.h"
#include "tsp/task.h"

#include <cstdint>
#include <iostream>
#include <memory>
#include <string>
#include <sys/types.h>


// w_ij --- distances
// c_ij --- is edge taken
// \sum_{ij} w_ij * c_ij -> min
// if i != j: 0 <= c_ij <= 1 (taken or not)
// if i == j: c_ij = 0
// \sum_{i} c_ij <= 1 (in-deg)
// \sum_{j} c_ij <= 1 (out-deg)
// 0 <= c_ij + c_ji <= 1 (only one direction)
void Solve(const tsp::AdjacencyMatrix& matrix) {
  using operations_research::MPConstraint;
  using operations_research::MPObjective;
  using operations_research::MPSolver;
  using operations_research::MPVariable;

  std::unique_ptr<MPSolver> solver(MPSolver::CreateSolver("SCIP"));
  if (!solver) {
    LOG(WARNING) << "SCIP solver unavailable.";
    return;
  }

  const std::vector<std::vector<std::string>> names = [&]() {
    std::vector<std::vector<std::string>> names(matrix.size(),
                                                std::vector<std::string>(matrix.size()));
    for (uint32_t i = 0; i < matrix.size(); ++i) {
      for (uint32_t j = 0; j < matrix.size(); ++j) {
        names[i][j] = "c_" + std::to_string(i) + "_" + std::to_string(j);
      }
    }
    return names;
  }();

  std::cerr << "OK" << std::endl;
  const std::vector<std::vector<MPVariable*>> is_edge_taken = [&]() {
    std::vector<std::vector<MPVariable*>> result(matrix.size(),
                                                 std::vector<MPVariable*>(matrix.size()));
    for (uint32_t i = 0; i < matrix.size(); ++i) {
      for (uint32_t j = 0; j < matrix.size(); ++j) {
        if (i != j) {
          // if i != j: 0 <= c_ij <= 1 (taken or not)
          result[i][j] = solver->MakeNumVar(0.0, 1.0, names[i][j]);
        } else {
          // if i == j: c_ij = 0
          result[i][j] = solver->MakeNumVar(0.0, 0.0, names[i][j]);
        }
      }
    }
    return result;
  }();
  LOG(INFO) << "Number of variables = " << solver->NumVariables();

  // 0 <= c_ij <= 1 (taken or not)
  for (uint32_t i = 0; i < matrix.size(); ++i) {
    for (uint32_t j = 0; j < matrix.size(); ++j) {
      MPConstraint* const c = solver->MakeRowConstraint(0, 1);
      c->SetCoefficient(is_edge_taken[i][j], 1);
    }
  }

  // \sum_{i} c_ij = 1 (in-deg)
  for (uint32_t j = 0; j < matrix.size(); ++j) {
    MPConstraint* const c = solver->MakeRowConstraint(1, 1);
    for (uint32_t i = 0; i < matrix.size(); ++i) {
      c->SetCoefficient(is_edge_taken[i][j], 1);
    }
  }

  // \sum_{j} c_ij = 1 (out-deg)
  for (uint32_t i = 0; i < matrix.size(); ++i) {
    MPConstraint* const c = solver->MakeRowConstraint(1, 1);
    for (uint32_t j = 0; j < matrix.size(); ++j) {
      c->SetCoefficient(is_edge_taken[i][j], 1);
    }
  }

  // 0 <= c_ij + c_ji <= 1 (only one direction)
  for (uint32_t i = 0; i < matrix.size(); ++i) {
    for (uint32_t j = i + 1; j < matrix.size(); ++j) {
      MPConstraint* const c = solver->MakeRowConstraint(0, 1);
      c->SetCoefficient(is_edge_taken[i][j], 1);
      c->SetCoefficient(is_edge_taken[j][i], 1);
    }
  }

  LOG(INFO) << "Number of constraints = " << solver->NumConstraints();

  // \sum_{ij} w_ij * c_ij -> min
  const MPObjective* objective = [&]() {
    MPObjective* const objective = solver->MutableObjective();
    for (uint32_t i = 0; i < matrix.size(); ++i) {
      for (uint32_t j = 0; j < matrix.size(); ++j) {
        objective->SetCoefficient(is_edge_taken[i][j], matrix[i][j]);
      }
    }
    objective->SetMinimization();
    return objective;
  }();

  const MPSolver::ResultStatus result_status = solver->Solve();
  // Check that the problem has an optimal solution.
  if (result_status != MPSolver::OPTIMAL) {
    LOG(FATAL) << "The problem does not have an optimal solution!";
  }

  LOG(INFO) << "Solution:";
  LOG(INFO) << "Optimal objective value = " << objective->Value();
}

ABSL_FLAG(std::string, input, {}, "Input file");

int main(int argc, char** argv) {
  absl::ParseCommandLine(argc, argv);
  std::string input = absl::GetFlag(FLAGS_input);

  tsp::AdjacencyMatrix matrix = tsp::AdjacencyMatrixFromFile(input);
  Solve(matrix);
  return 0;
}
