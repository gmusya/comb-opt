#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace tsp {

using Weight = uint32_t;
using AdjacencyMatrix = std::vector<std::vector<Weight>>;
using Vertex = uint32_t;
using Solution = std::vector<Vertex>;

AdjacencyMatrix AdjacencyMatrixFromFile(const std::string& filename);

Solution SolutionFromFile(const std::string& filename);

void SolutionToFile(const Solution& solution, const std::string& filename);

Weight GetScore(const AdjacencyMatrix& matrix, const Solution& solution);

}// namespace tsp
