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

template<typename Callback>
void DFS(int v, std::vector<bool>& used, const tsp::AdjacencyLists& graph, Callback& callback) {
  used[v] = true;
  callback(v);
  for (uint32_t u : graph[v]) {
    if (used[u]) {
      continue;
    }
    DFS(u, used, graph, callback);
  }
}

// w_ij --- distances
// c_ij --- is edge taken
// \sum_{ij} w_ij * c_ij -> min
// if i != j: 0 <= c_ij <= 1 (taken or not)
// if i == j: c_ij = 0
// \sum_{i} c_ij <= 1 (in-deg)
// \sum_{j} c_ij <= 1 (out-deg)
// 0 <= c_ij + c_ji <= 1 (only one direction)
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

  while (true) {
    const MPSolver::ResultStatus result_status = solver->Solve();

    // Check that the problem has an optimal solution.
    if (result_status != MPSolver::OPTIMAL) {
      LOG(FATAL) << "The problem does not have an optimal solution!";
    }

    LOG(INFO) << "Solution:";
    LOG(INFO) << "Optimal objective value = " << objective->Value();

    tsp::AdjacencyLists adj_lists(matrix.size());
    for (uint32_t i = 0; i < matrix.size(); ++i) {
      for (uint32_t j = 0; j < matrix.size(); ++j) {
        if (is_edge_taken[i][j]->solution_value() == 1) {
          adj_lists[i].push_back(j);
          adj_lists[j].push_back(i);
        }
      }
    }

    std::vector<std::vector<tsp::Vertex>> comps;
    std::vector<bool> used(matrix.size());
    for (tsp::Vertex v = 0; v < matrix.size(); ++v) {
      if (!used[v]) {
        comps.emplace_back();
        auto lambda = [&](tsp::Vertex v) {
          comps.back().emplace_back(v);
        };
        DFS(v, used, adj_lists, lambda);
      }
    }

    if (comps.size() == 1) {
      return comps[0];
    }

    for (const auto& comp : comps) {
      MPConstraint* const c = solver->MakeRowConstraint(0, comp.size() - 1);
      for (tsp::Vertex v : comp) {
        for (tsp::Vertex u : comp) {
          c->SetCoefficient(is_edge_taken[v][u], 1);
        }
      }
    }
  }
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
