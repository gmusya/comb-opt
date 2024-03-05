#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "tsp/task.h"

#include <cstdint>
#include <google/protobuf/stubs/port.h>
#include <limits>
#include <random>
#include <vector>

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

std::vector<tsp::Vertex> FindOptimalPath(const tsp::AdjacencyMatrix& matrix,
                                         const std::vector<tsp::Vertex>& vertices) {
  uint32_t size = vertices.size();
  constexpr tsp::Weight kInfinity = std::numeric_limits<tsp::Weight>::max();
  std::vector<std::vector<tsp::Weight>> result(size,
                                               std::vector<tsp::Weight>(1 << size, kInfinity));
  std::vector<std::vector<tsp::Vertex>> prev(size, std::vector<tsp::Vertex>(1 << size, 0));
  std::vector<uint32_t> masks(1 << size);
  std::iota(masks.begin(), masks.end(), 0);
  std::sort(masks.begin(), masks.end(), [](uint32_t lhs, uint32_t rhs) {
    return __builtin_popcount(lhs) < __builtin_popcount(rhs);
  });

  result[0][1] = 0;
  for (uint32_t mask : masks) {
    for (uint32_t last = 0; last < size; ++last) {
      if (result[last][mask] == kInfinity) {
        continue;
      }
      for (uint32_t after = 0; after < size; ++after) {
        uint32_t new_mask = mask | (1u << after);
        if (mask == new_mask) {
          continue;
        }
        if (result[last][mask] + matrix[vertices[last]][vertices[after]] <
            result[after][new_mask]) {
          result[after][new_mask] = result[last][mask] + matrix[vertices[last]][vertices[after]];
          prev[after][new_mask] = last;
        }
      }
    }
  }

  std::vector<tsp::Vertex> optimal_path;
  tsp::Vertex last = size - 1;
  uint32_t last_mask = (1u << size) - 1;
  // std::cerr << "Value = " << result[last][last_mask] << std::endl;
  while (last != 0) {
    optimal_path.emplace_back(vertices[last]);
    auto new_mask = (last_mask ^ (1u << last));
    last = prev[last][last_mask];
    last_mask = new_mask;
  }
  optimal_path.emplace_back(vertices[last]);

  return optimal_path;
}

std::pair<tsp::Solution, bool> ImproveNeighbours(const tsp::AdjacencyMatrix& matrix,
                                                 const tsp::Solution initial_solution,
                                                 uint32_t max_permutation_size) {
  auto solution = initial_solution;
  tsp::Weight solution_score = tsp::GetScore(matrix, solution);
  bool is_improved = false;
  for (uint32_t permutation_size = 2; permutation_size <= max_permutation_size;
       ++permutation_size) {
    // std::cerr << "ImproveNeighbours, permutation_size = " << permutation_size << std::endl;

    for (uint32_t i = 0; i < matrix.size(); ++i) {
      // std::cerr << "ImproveNeighbours, permutation_size = " << permutation_size << "(" << i << "/"
      // << matrix.size() << ")" << std::endl;

      std::vector<tsp::Vertex> vertices;
      for (uint32_t j = 0; j < permutation_size; ++j) {
        vertices.emplace_back(solution[(i + j) % matrix.size()]);
      }
      vertices = FindOptimalPath(matrix, vertices);

      auto new_solution = solution;
      for (uint32_t j = 0; j < permutation_size; ++j) {
        new_solution[(i + j) % matrix.size()] = vertices[j];
      }
      auto new_score = tsp::GetScore(matrix, new_solution);
      if (new_score < solution_score) {
        solution_score = new_score;
        solution = new_solution;
        is_improved = true;
        // std::cerr << "score = " << solution_score << std::endl;
        break;
      }
    }
    if (is_improved) {
      return std::make_pair(solution, is_improved);
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
    // std::cerr << "ImproveAny, permutation_size = " << permutation_size << std::endl;

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
  return std::make_pair(solution, is_improved);
}

std::pair<tsp::Solution, bool> ImproveShift(const tsp::AdjacencyMatrix& matrix,
                                            const tsp::Solution initial_solution) {
  auto solution = initial_solution;
  tsp::Weight solution_score = tsp::GetScore(matrix, solution);
  bool is_improved = false;
  std::cerr << "ImproveShift" << std::endl;
  for (uint32_t from = 0; from < matrix.size(); ++from) {
    std::cerr << "ImproveShift (" << from << "/" << matrix.size() << ")" << std::endl;
    for (uint32_t len = 1; len < matrix.size() / 2; ++len) {
      std::vector<tsp::Vertex> old_v;
      std::vector<tsp::Vertex> new_v;
      for (uint32_t sh = 0; sh < matrix.size(); ++sh) {
        uint32_t ind = (from + sh) % matrix.size();
        if (sh < len) {
          old_v.emplace_back(solution[ind]);
        } else {
          new_v.emplace_back(solution[ind]);
        }
      }

      for (uint32_t to = 0; to <= new_v.size(); ++to) {
        std::vector<tsp::Vertex> new_solution = new_v;
        new_solution.insert(new_solution.begin() + to, old_v.begin(), old_v.end());

        auto [new_new_solution, is_improved_2] = ImproveNeighbours(matrix, new_solution, 5);
        tsp::Weight new_score = tsp::GetScore(matrix, new_new_solution);
        new_solution = std::move(new_new_solution);

        if (new_score < solution_score) {
          solution_score = new_score;
          solution = new_solution;
          is_improved = true;
          std::cerr << "new score = " << tsp::GetScore(matrix, solution) << std::endl;
          return std::make_pair(solution, is_improved);
        }
      }
    }
  }
  return std::make_pair(solution, is_improved);
}

tsp::Solution LocalSolve(const tsp::AdjacencyMatrix& matrix,
                         const tsp::Solution& initial_solution) {
  auto solution = initial_solution;
  while (true) {
    {
      auto [new_solution, is_improved] = ImproveNeighbours(matrix, solution, 9);
      if (is_improved) {
        solution = new_solution;
        continue;
      }
    }
    {
      auto [new_solution, is_improved] = ImproveShift(matrix, solution);
      if (is_improved) {
        solution = new_solution;
        continue;
      }
    }
    {
      auto [new_solution, is_improved] = ImproveAny(matrix, solution, 4);
      if (is_improved) {
        solution = new_solution;
        continue;
      }
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
