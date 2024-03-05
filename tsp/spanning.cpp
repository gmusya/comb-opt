#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "tsp/task.h"
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

using Edge = std::pair<tsp::Vertex, tsp::Vertex>;

class Dsu {
  public:
  explicit Dsu(uint32_t vertices)
      : representative_([](uint32_t vertices) {
          std::vector<tsp::Vertex> result(vertices);
          std::iota(result.begin(), result.end(), 0);
          return result;
        }(vertices)),
        rank_(vertices) {
  }

  tsp::Vertex GetRepresentative(tsp::Vertex v) const {
    if (v == representative_[v]) {
      return v;
    } else {
      return representative_[v] = GetRepresentative(representative_[v]);
    }
  }

  bool Merge(tsp::Vertex u, tsp::Vertex v) {
    u = GetRepresentative(u);
    v = GetRepresentative(v);
    if (u == v) {
      return false;
    }
    if (rank_[u] > rank_[v]) {
      std::swap(u, v);
    }
    if (rank_[u] == rank_[v]) {
      ++rank_[v];
    }
    representative_[u] = v;
    return true;
  }

  private:
  mutable std::vector<tsp::Vertex> representative_;
  std::vector<uint32_t> rank_;
};

std::vector<Edge> MinimalSpanningTree(const tsp::AdjacencyMatrix& graph) {
  std::vector<Edge> all_edges;
  all_edges.reserve(graph.size() * (graph.size() - 1) / 2);
  for (tsp::Vertex i = 0; i < graph.size(); ++i) {
    for (tsp::Vertex j = i + 1; j < graph.size(); ++j) {
      all_edges.emplace_back(i, j);
    }
  }
  std::sort(all_edges.begin(), all_edges.end(), [&graph](const auto& lhs, const auto& rhs) {
    return graph[lhs.first][lhs.second] < graph[rhs.first][rhs.second];
  });
  std::vector<Edge> spanning_tree;
  Dsu dsu(graph.size());
  for (const auto& [u, v] : all_edges) {
    if (dsu.Merge(u, v)) {
      spanning_tree.emplace_back(u, v);
    }
  }
  if (spanning_tree.size() + 1 != graph.size()) {
    throw std::runtime_error("Graph is not connected");
  }
  return spanning_tree;
}

tsp::AdjacencyLists GraphFromEdges(const std::vector<Edge>& edges) {
  uint32_t max_value = 0;
  for (const auto& [u, v] : edges) {
    max_value = std::max({max_value, u, v});
  }
  tsp::AdjacencyLists g(max_value + 1);
  for (const auto& [u, v] : edges) {
    g[u].emplace_back(v);
    g[v].emplace_back(u);
  }
  return g;
}

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

std::vector<tsp::Vertex> EulerFromGraph(const tsp::AdjacencyLists& graph) {
  std::vector<bool> used(graph.size());
  std::vector<tsp::Vertex> result;
  auto callback = [&result](uint32_t v) mutable {
    result.emplace_back(v);
  };
  DFS(0, used, graph, callback);
  if (result.size() != graph.size()) {
    throw std::runtime_error("Graph is not connected");
  }
  return result;
}

tsp::Solution SpanningSolve(const tsp::AdjacencyMatrix& matrix) {
  std::vector<Edge> spanning_tree = MinimalSpanningTree(matrix);
  uint32_t spanning_tree_path_size = 0;
  for (const auto& [u, v] : spanning_tree) {
    spanning_tree_path_size += matrix[u][v];
  }
  std::cerr << "spanning tree double weight = " << 2 * spanning_tree_path_size << std::endl;
  std::vector<std::vector<uint32_t>> tree = GraphFromEdges(spanning_tree);
  return EulerFromGraph(tree);
}

ABSL_FLAG(std::string, input, {}, "Input file");
ABSL_FLAG(std::string, output, {}, "Output file");

int main(int argc, char** argv) {
  absl::ParseCommandLine(argc, argv);
  std::string input = absl::GetFlag(FLAGS_input);
  std::string output = absl::GetFlag(FLAGS_output);
  tsp::AdjacencyMatrix matrix = tsp::AdjacencyMatrixFromFile(input);

  tsp::Solution solution = SpanningSolve(matrix);

  tsp::Weight score = tsp::GetScore(matrix, solution);
  std::cerr << "score = " << score << std::endl;
  tsp::SolutionToFile(solution, output);
  return 0;
}
