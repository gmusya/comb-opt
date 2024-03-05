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
ABSL_FLAG(std::string, initial_route, {}, "File with route to optimize");
ABSL_FLAG(uint32_t, seed, 2101, "Seed for random");

std::pair<tsp::Solution, bool> ImproveNeighbours(const tsp::AdjacencyMatrix& matrix,
                                                 const tsp::Solution initial_solution,
                                                 uint32_t max_permutation_size) {
  auto solution = initial_solution;
  tsp::Weight solution_score = tsp::GetScore(matrix, solution);
  bool is_improved = false;
  for (uint32_t permutation_size = 2; permutation_size <= max_permutation_size;
       ++permutation_size) {
    std::cerr << "ImproveNeighbours, permutation_size = " << permutation_size << std::endl;

    for (uint32_t i = 0; i < matrix.size(); ++i) {
      std::cerr << "ImproveNeighbours, permutation_size = " << permutation_size << "(" << i << "/"
                << matrix.size() << ")" << std::endl;

      std::vector<uint32_t> perm(permutation_size);
      std::iota(perm.begin(), perm.end(), 0);
      do {
        auto new_solution = solution;
        for (uint32_t j = 0; j < permutation_size; ++j) {
          new_solution[(i + j) % matrix.size()] = solution[(i + perm[j]) % matrix.size()];
        }
        tsp::Weight new_score = tsp::GetScore(matrix, new_solution);
        if (new_score < solution_score) {
          solution_score = new_score;
          solution = new_solution;
          is_improved = true;
          std::cerr << "new score = " << tsp::GetScore(matrix, solution) << std::endl;
          break;
        }
      } while (std::next_permutation(perm.begin(), perm.end()));
    }
  }
  return std::make_pair(solution, is_improved);
}

std::pair<tsp::Solution, bool> ImproveAny(const tsp::AdjacencyMatrix& matrix,
                                          const tsp::Solution initial_solution,
                                          uint32_t max_permutation_size) {
  auto solution = initial_solution;
  tsp::Weight solution_score = tsp::GetScore(matrix, solution);
  bool is_improved = false;
  for (uint32_t permutation_size = 2; permutation_size <= max_permutation_size;
       ++permutation_size) {
    std::cerr << "ImproveAny, permutation_size = " << permutation_size << std::endl;

    for (uint32_t i = 0; i < matrix.size(); ++i) {
      std::cerr << "ImproveAny, permutation_size = " << permutation_size << "(" << i << "/"
                << matrix.size() << ")" << std::endl;

      std::vector<uint32_t> to_permute(matrix.size(), 0);
      for (uint32_t i = matrix.size() - 1; i >= matrix.size() - permutation_size; --i) {
        to_permute[i] = 1;
      }
      do {
        std::vector<uint32_t> ones_positions;
        for (uint32_t i = 0; i < matrix.size(); ++i) {
          if (to_permute[i] == 1) {
            ones_positions.push_back(i);
          }
        }
        std::vector<uint32_t> perm(permutation_size);
        std::iota(perm.begin(), perm.end(), 0);
        do {
          auto new_solution = solution;
          for (uint32_t j = 0; j < permutation_size; ++j) {
            new_solution[ones_positions[j]] = solution[ones_positions[perm[j]]];
          }

          tsp::Weight new_score = tsp::GetScore(matrix, new_solution);
          if (new_score < solution_score) {
            solution_score = new_score;
            solution = new_solution;
            is_improved = true;
            std::cerr << "new score = " << tsp::GetScore(matrix, solution) << std::endl;
            break;
          }
        } while (std::next_permutation(perm.begin(), perm.end()));
      } while (std::next_permutation(to_permute.begin(), to_permute.end()));
    }
  }
  return std::make_pair(solution, is_improved);
}


tsp::Solution LocalSolve(const tsp::AdjacencyMatrix& matrix,
                         const tsp::Solution& initial_solution) {
  auto solution = initial_solution;
  while (true) {
    {
      auto [new_solution, is_improved] = ImproveNeighbours(matrix, solution, 8);
      if (is_improved) {
        solution = new_solution;
        continue;
      }
    }
    auto [new_solution, is_improved] = ImproveAny(matrix, solution, 3);
    if (is_improved) {
      solution = new_solution;
      continue;
    }
    break;
  }
  return solution;
}

int main(int argc, char** argv) {
  absl::ParseCommandLine(argc, argv);
  std::string input = absl::GetFlag(FLAGS_input);
  std::string output = absl::GetFlag(FLAGS_output);
  std::string initial_route = absl::GetFlag(FLAGS_initial_route);
  uint32_t seed = absl::GetFlag(FLAGS_seed);

  tsp::AdjacencyMatrix matrix = tsp::AdjacencyMatrixFromFile(input);
  auto solution = tsp::SolutionFromFile(initial_route);
  std::cerr << "initial score = " << tsp::GetScore(matrix, solution) << std::endl;

  solution = LocalSolve(matrix, solution);
  tsp::Weight score = tsp::GetScore(matrix, solution);
  std::cerr << "score = " << score << std::endl;
  tsp::SolutionToFile(solution, output);
  return 0;
}
