#include "absl/flags/flag.h"
#include "absl/flags/parse.h"

#include <iostream>
#include <string>

ABSL_FLAG(std::string, input, {}, "Input file");
ABSL_FLAG(std::string, output, {}, "Output file");

int main(int argc, char **argv) {
  absl::ParseCommandLine(argc, argv);
  std::string input = absl::GetFlag(FLAGS_input);
  std::string output = absl::GetFlag(FLAGS_output);
  std::cout << input << ' ' << output << std::endl;
  return 0;
}
