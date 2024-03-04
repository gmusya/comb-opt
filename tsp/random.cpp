#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "tsp/task.h"

#include <cstdint>
#include <random>

tsp::Solution RandomSolve(const tsp::AdjacencyMatrix& matrix, uint32_t seed) {
  uint32_t vertices_count = matrix.size();
  tsp::Solution result(vertices_count);
  std::iota(result.begin(), result.end(), 0);
  std::mt19937 rnd(seed);
  std::shuffle(result.begin(), result.end(), rnd);
  return result;
}

ABSL_FLAG(std::string, input, {}, "Input file");
ABSL_FLAG(std::string, output, {}, "Output file");
ABSL_FLAG(uint32_t, seed, 2101, "Seed for random");

int main(int argc, char** argv) {
  absl::ParseCommandLine(argc, argv);
  std::string input = absl::GetFlag(FLAGS_input);
  std::string output = absl::GetFlag(FLAGS_output);
  uint32_t seed = absl::GetFlag(FLAGS_seed);

  tsp::AdjacencyMatrix matrix = tsp::AdjacencyMatrixFromFile(input);
  tsp::Solution solution = RandomSolve(matrix, seed);
  tsp::Weight score = tsp::GetScore(matrix, solution);
  std::cerr << "score = " << score << std::endl;
  tsp::SolutionToFile(solution, output);
  return 0;
}
