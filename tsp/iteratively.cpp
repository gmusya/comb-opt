#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "tsp/task.h"

#include <cstdint>
#include <vector>

tsp::Solution IterativelySolve(const tsp::AdjacencyMatrix& matrix) {
  std::vector<tsp::Vertex> solution = {0, 1};

  for (tsp::Vertex v = 2; v < matrix.size(); ++v) {
    std::vector<std::vector<tsp::Vertex>> new_solutions = {};
    for (uint32_t pos = 0; pos < solution.size(); ++pos) {
      auto new_solution = solution;
      new_solution.insert(new_solution.begin() + pos, v);
      new_solutions.emplace_back(std::move(new_solution));
    }

    std::sort(new_solutions.begin(), new_solutions.end(),
              [&matrix](const auto& lhs, const auto& rhs) {
                return tsp::GetScore(matrix, lhs) < tsp::GetScore(matrix, rhs);
              });

    solution = new_solutions[0];
  }

  return solution;
}

ABSL_FLAG(std::string, input, {}, "Input file");
ABSL_FLAG(std::string, output, {}, "Output file");

int main(int argc, char** argv) {
  absl::ParseCommandLine(argc, argv);
  std::string input = absl::GetFlag(FLAGS_input);
  std::string output = absl::GetFlag(FLAGS_output);

  tsp::AdjacencyMatrix matrix = tsp::AdjacencyMatrixFromFile(input);

  auto solution = IterativelySolve(matrix);

  tsp::Weight score = tsp::GetScore(matrix, solution);
  std::cerr << "score = " << score << std::endl;
  tsp::SolutionToFile(solution, output);
  return 0;
}
