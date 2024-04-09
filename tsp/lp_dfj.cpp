#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "ortools/linear_solver/linear_solver.h"
#include "tsp/task.h"

#include <cstdint>
#include <fstream>
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

using operations_research::MPConstraint;
using operations_research::MPObjective;
using operations_research::MPSolver;
using operations_research::MPVariable;

void PrintState(const std::vector<std::vector<MPVariable*>>& is_edge_taken,
                std::string output_file) {
  static int iter = 0;
  ++iter;
  std::ofstream out(output_file + std::to_string(iter) + ".txt");
  out << is_edge_taken.size() << '\n';
  for (int i = 0; i < is_edge_taken.size(); ++i) {
    for (int j = 0; j < is_edge_taken.size(); ++j) {
      out << is_edge_taken[i][j]->solution_value() << ' ';
    }
    out << '\n';
  }
}

tsp::AdjacencyLists GetBlackEdges(const std::vector<std::vector<MPVariable*>>& is_edge_taken) {
  tsp::AdjacencyLists adj_lists(is_edge_taken.size());
  for (uint32_t i = 0; i < is_edge_taken.size(); ++i) {
    for (uint32_t j = 0; j < is_edge_taken.size(); ++j) {
      if (is_edge_taken[i][j]->solution_value() == 1) {
        adj_lists[i].push_back(j);
        adj_lists[j].push_back(i);
      }
    }
  }
  return adj_lists;
}

tsp::AdjacencyLists GetRedEdges(const std::vector<std::vector<MPVariable*>>& is_edge_taken) {
  tsp::AdjacencyLists adj_lists(is_edge_taken.size());
  for (uint32_t i = 0; i < is_edge_taken.size(); ++i) {
    for (uint32_t j = 0; j < is_edge_taken.size(); ++j) {
      if (is_edge_taken[i][j]->solution_value() == 0.5) {
        adj_lists[i].push_back(j);
        adj_lists[j].push_back(i);
      }
    }
  }
  return adj_lists;
}


std::vector<std::vector<tsp::Vertex>> GetComps(const tsp::AdjacencyLists& adj_lists) {
  std::vector<std::vector<tsp::Vertex>> comps;
  std::vector<bool> used(adj_lists.size());
  for (tsp::Vertex v = 0; v < adj_lists.size(); ++v) {
    if (!used[v]) {
      comps.emplace_back();
      auto lambda = [&](tsp::Vertex v) {
        comps.back().emplace_back(v);
      };
      DFS(v, used, adj_lists, lambda);
    }
  }
  return comps;
}

void RemoveBlackCycles(MPSolver* solver, const MPObjective* objective,
                       const std::vector<std::vector<MPVariable*>>& is_edge_taken,
                       const std::string& output) {
  static int calls = 0;
  static int cycles = 0;
  const MPSolver::ResultStatus result_status = solver->Solve();
  if (result_status != MPSolver::OPTIMAL) {
    LOG(FATAL) << "The problem does not have an optimal solution!";
  }
  std::cerr << "RemoveBlackCycles start: " << objective->Value() << std::endl;
  tsp::AdjacencyLists old_lists;

  while (true) {
    std::cerr << "Solve start" << std::endl;
    const MPSolver::ResultStatus result_status = solver->Solve();
    if (result_status != MPSolver::OPTIMAL) {
      LOG(FATAL) << "The problem does not have an optimal solution!";
    }
    std::cerr << "Score: " << objective->Value() << std::endl;
    PrintState(is_edge_taken, output);

    auto adj_lists = GetBlackEdges(is_edge_taken);
    if (old_lists == adj_lists) {
      break;
    }
    old_lists = adj_lists;

    std::vector<std::vector<tsp::Vertex>> comps = GetComps(adj_lists);
    if (comps.size() == 1) {
      break;
    }
    ++calls;
    for (const auto& comp : comps) {
      MPConstraint* const c = solver->MakeRowConstraint(0, comp.size() - 1);
      ++cycles;
      for (tsp::Vertex v : comp) {
        for (tsp::Vertex u : comp) {
          c->SetCoefficient(is_edge_taken[v][u], 1);
        }
      }
    }
  }

  std::cerr << "calls = " << calls << ", cycles = " << cycles << std::endl;
}

std::vector<std::vector<tsp::Vertex>>
GetTeeth(const std::vector<std::vector<MPVariable*>>& is_edge_taken) {
  std::vector<std::vector<tsp::Vertex>> teeth;
  auto black_edges = GetBlackEdges(is_edge_taken);
  return GetComps(black_edges);
}

std::vector<std::vector<tsp::Vertex>>
GetPossibleHandles(const std::vector<std::vector<MPVariable*>>& is_edge_taken) {
  std::vector<std::vector<tsp::Vertex>> teeth;
  auto black_edges = GetRedEdges(is_edge_taken);
  return GetComps(black_edges);
}

// w_ij --- distances
// c_ij --- is edge taken
// \sum_{ij} w_ij * c_ij -> min
// if i < j: 0 <= c_ij <= 1 (taken or not)
// if i == j: c_ij = 0
// \sum_{i} c_ij <= 2
void LPLowerBound(const tsp::AdjacencyMatrix& matrix, std::string output) {
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

  // if i < j: 0 <= c_ij <= 1 (taken or not)
  for (uint32_t i = 0; i < matrix.size(); ++i) {
    for (uint32_t j = i + 1; j < matrix.size(); ++j) {
      MPConstraint* const c = solver->MakeRowConstraint(0, 1);
      c->SetCoefficient(is_edge_taken[i][j], 1);
    }
  }
  // \sum_{i} c_ij = 2
  for (uint32_t j = 0; j < matrix.size(); ++j) {
    MPConstraint* const c = solver->MakeRowConstraint(2, 2);
    for (uint32_t i = 0; i < matrix.size(); ++i) {
      if (i == j) {
        continue;
      }
      int x = std::min(i, j);
      int y = std::max(i, j);
      c->SetCoefficient(is_edge_taken[x][y], 1);
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


  int cnt = 0;
  output = output.substr(0, output.size() - 4);

  // Try to find violated combs
  while (true) {
    RemoveBlackCycles(solver.get(), objective, is_edge_taken, output);
    auto teeth = GetTeeth(is_edge_taken);
    auto maybe_handles = GetPossibleHandles(is_edge_taken);

    std::vector<int32_t> from_vertex_to_tooth(matrix.size(), -1);
    for (size_t i = 0; i < teeth.size(); ++i) {
      for (const auto& v : teeth[i]) {
        from_vertex_to_tooth[v] = i;
      }
    }

    std::vector<std::vector<double>> old_solution_values(matrix.size(),
                                                         std::vector<double>(matrix.size()));
    for (tsp::Vertex x = 0; x < matrix.size(); ++x) {
      for (tsp::Vertex y = 0; y < matrix.size(); ++y) {
        old_solution_values[x][y] = is_edge_taken[x][y]->solution_value();
      }
    }

    int changes = 0;
    static int violated_combs = 0;
    for (const auto& handle : maybe_handles) {
      if (handle.size() % 2 == 0 || handle.size() == 1) {
        continue;
      }
      ++changes;
      std::set<tsp::Vertex> handle_vertices(handle.begin(), handle.end());

      std::set<int32_t> tooth_ids;
      for (tsp::Vertex v : handle) {
        if (from_vertex_to_tooth[v] != -1) {
          int32_t id = from_vertex_to_tooth[v];
          tooth_ids.insert(id);
        }
      }

      std::map<std::pair<tsp::Vertex, tsp::Vertex>, int32_t> edges_to_ban;
      for (tsp::Vertex v : handle) {
        for (tsp::Vertex u = 0; u < matrix.size(); ++u) {
          if (handle_vertices.contains(u)) {
            continue;
          }
          edges_to_ban[std::make_pair(std::min(u, v), std::max(u, v))]++;
        }
      }

      for (int32_t id : tooth_ids) {
        for (tsp::Vertex v : teeth[id]) {
          for (tsp::Vertex u = 0; u < matrix.size(); ++u) {
            if (from_vertex_to_tooth[u] == id) {
              continue;
            }
            edges_to_ban[std::make_pair(std::min(u, v), std::max(u, v))]++;
          }
        }
      }

      int32_t bound = 3 * tooth_ids.size() + 1;
      std::cerr << "bound = " << bound << std::endl;
      ++violated_combs;
      std::cerr << "violated_combs = " << violated_combs << std::endl;
      MPConstraint* const qwe = solver->MakeRowConstraint(bound, solver->infinity());
      for (const auto& [key, cnt] : edges_to_ban) {
        const auto& [x, y] = key;
        qwe->SetCoefficient(is_edge_taken[x][y], cnt);
      }
    }

    const MPSolver::ResultStatus result_status = solver->Solve();

    if (result_status != MPSolver::OPTIMAL) {
      LOG(FATAL) << "The problem does not have an optimal solution!";
    }

    std::cerr << "LB: " << objective->Value() << std::endl;
    if (changes == 0) {
      break;
    }
  }
  std::cerr << "Result: " << objective->Value() << std::endl;
}

ABSL_FLAG(std::string, input, {}, "Input file");
ABSL_FLAG(std::string, output, {}, "Output file");

int main(int argc, char** argv) {
  absl::ParseCommandLine(argc, argv);
  std::string input = absl::GetFlag(FLAGS_input);
  std::string output = absl::GetFlag(FLAGS_output);

  tsp::AdjacencyMatrix matrix = tsp::AdjacencyMatrixFromFile(input);

  LPLowerBound(matrix, output);

  return 0;
}
