#include "tsp/task.h"

#include <cassert>
#include <cstdint>
#include <fstream>

namespace tsp {

AdjacencyMatrix AdjacencyMatrixFromFile(const std::string& filename) {
  std::ifstream input(filename);
  uint32_t vertices_count;
  input >> vertices_count;
  AdjacencyMatrix result(vertices_count, std::vector<uint32_t>(vertices_count));
  for (Vertex u = 0; u < vertices_count; ++u) {
    for (Vertex v = 0; v < vertices_count; ++v) {
      input >> result[u][v];
    }
  }
  return result;
}

Solution SolutionFromFile(const std::string& filename) {
  std::ifstream input(filename);
  uint32_t vertices_count;
  input >> vertices_count;
  Solution result(vertices_count);
  for (uint32_t i = 0; i < vertices_count; ++i) {
    input >> result[i];
  }
  return result;
}

void SolutionToFile(const Solution& solution, const std::string& filename) {
  std::ofstream output(filename);
  uint32_t vertices_count = solution.size();
  output << vertices_count << '\n';
  for (uint32_t i = 0; i < vertices_count; ++i) {
    if (i != 0) {
      output << ' ';
    }
    output << solution[i];
  }
  output << '\n';
}

Weight GetScore(const AdjacencyMatrix& matrix, const Solution& solution) {
  assert(matrix.size() == solution.size());
  Weight result = 0;
  for (uint32_t i = 0; i < solution.size(); ++i) {
    Vertex u = solution[i];
    Vertex v = (i + 1 == solution.size()) ? solution[0] : solution[i + 1];
    result += matrix[u][v];
  }
  return result;
}

}// namespace tsp
