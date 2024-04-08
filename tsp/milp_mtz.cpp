#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "ortools/linear_solver/linear_solver.h"
#include "tsp/task.h"

#include <cstdint>
#include <iostream>
#include <memory>
#include <stdexcept>
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
// o_i --- position of i-th vertex in order
// o_i - o_j + 1 <= (n - 1)(1 - c_ij) (if c_ij = 1 then o_i < o_j)
tsp::Solution MILPSolve(const tsp::AdjacencyMatrix& matrix) {
  using operations_research::MPConstraint;
  using operations_research::MPObjective;
  using operations_research::MPSolver;
  using operations_research::MPVariable;

  std::unique_ptr<MPSolver> solver(MPSolver::CreateSolver("SCIP"));
  if (!solver) {
    LOG(ERROR) << "SCIP solver unavailable.";
    throw std::runtime_error("SCIP solver unavailable.");
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
          result[i][j] = solver->MakeIntVar(0.0, 1.0, names[i][j]);
        } else {
          // if i == j: c_ij = 0
          result[i][j] = solver->MakeIntVar(0.0, 0.0, names[i][j]);
        }
      }
    }
    return result;
  }();

  const std::vector<MPVariable*> order = [&]() {
    std::vector<MPVariable*> result(matrix.size());
    for (uint32_t i = 0; i < matrix.size(); ++i) {
      if (i == 0) {
        result[i] = solver->MakeIntVar(1, 1, "o_1");
      } else {
        result[i] = solver->MakeIntVar(2, matrix.size(), "o_" + std::to_string(i));
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

  // o_i - o_j + 1 <= (n - 1)(1 - c_ij) (if c_ij = 1 then o_i < o_j)
  // o_i - o_j + 1 <= (n - 1) - (n - 1)(c_ij)
  // o_i - o_j + (n - 1) c_ij <= n - 2
  for (uint32_t i = 1; i < matrix.size(); ++i) {
    for (uint32_t j = 1; j < matrix.size(); ++j) {
      if (i == j) {
        continue;
      }
      MPConstraint* const c =
              solver->MakeRowConstraint(-static_cast<int32_t>(matrix.size()), matrix.size() - 2);
      c->SetCoefficient(order[i], 1);
      c->SetCoefficient(order[j], -1);
      c->SetCoefficient(is_edge_taken[i][j], matrix.size() - 1);
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

  tsp::Solution solution;
  uint32_t v = 0;
  do {
    uint32_t u = 0;
    while (u < matrix.size() && is_edge_taken[v][u]->solution_value() != 1) {
      ++u;
    }
    if (u == matrix.size()) {
      throw std::runtime_error("Internal error: incorrect solution");
    }
    solution.push_back(u);
    v = u;
  } while (v != 0);

  return solution;
}

ABSL_FLAG(std::string, input, {}, "Input file");
ABSL_FLAG(std::string, output, {}, "Output file");

int main(int argc, char** argv) {
  absl::ParseCommandLine(argc, argv);
  std::string input = absl::GetFlag(FLAGS_input);
  std::string output = absl::GetFlag(FLAGS_output);

  tsp::AdjacencyMatrix matrix = tsp::AdjacencyMatrixFromFile(input);

  tsp::Solution solution = MILPSolve(matrix);

  tsp::Weight score = tsp::GetScore(matrix, solution);
  std::cerr << "score = " << score << std::endl;
  tsp::SolutionToFile(solution, output, score);

  return 0;
}
