//
// Created by kellerberrin on 4/05/18.
//

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Setup logger and read the XML program options.

#include "kgl_gene_app.h"
#include "kel_utility.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>

// Define namespace alias
namespace kgl = kellerberrin::genome;
namespace kel = kellerberrin;


/// Prints a fatal error message to stderr and exits the program.
[[noreturn]] static void fatalExit(const std::string& message) {

  std::cerr << message << std::endl;
  std::cerr << kgl::GeneExecEnv::MODULE_NAME << " exits" << std::endl;
  std::exit(EXIT_FAILURE);

}


// Simple command-line parser replacing boost::program_options.
// Supports --flag=value and --flag value syntax for three required string flags
// plus --help/-h as a boolean flag.
struct ParsedArgs {

  std::string workDirectory;
  std::string logFile;
  std::string optionFile;
  bool helpRequested{false};

};


/// Parses command-line arguments into a ParsedArgs structure, supporting --flag=value and --flag value syntax.
[[nodiscard]] static std::optional<ParsedArgs> parseArgs(int argc, char const** argv) {

  ParsedArgs args;

  // Map flag names to the ParsedArgs members they populate.
  static const std::unordered_map<std::string, std::string*> string_flags = {
    {"--workDirectory", &args.workDirectory},
    {"--logFile",       &args.logFile},
    {"--optionFile",    &args.optionFile},
  };

  for (int i = 1; i < argc; ++i) {

    std::string arg = argv[i];

    if (arg == "--help" or arg == "-h") {

      args.helpRequested = true;
      continue;

    }

    // Handle --flag=value syntax.
    auto eq_pos = arg.find('=');
    std::string flag;
    std::string value;

    if (eq_pos != std::string::npos) {

      flag = arg.substr(0, eq_pos);
      value = arg.substr(eq_pos + 1);

    } else {

      // Handle --flag value syntax (next argument is the value).
      flag = arg;
      if (i + 1 >= argc) {

        std::cerr << "ERROR: missing value for " << flag << std::endl;
        return std::nullopt;

      }
      value = argv[++i];

    }

    auto it = string_flags.find(flag);
    if (it == string_flags.end()) {

      std::cerr << "ERROR: unknown option '" << flag << "'" << std::endl;
      return std::nullopt;

    }
    *it->second = value;

  }

  return args;

}


/// Parses the command line arguments and initializes the runtime environment.
bool kgl::GeneExecEnv::parseCommandLine(int argc, char const ** argv)
{

  std::stringstream ss;
  ss << "Population Genome Comparison, module: "
     << MODULE_NAME
     << " version: "
     << VERSION << '\n'
     << "Usage: --workDirectory=<work_directory> --logFile=<log_file.log> --optionFile=<option_file.xml> (all arguments required)";
  const std::string help_description = ss.str();

  if (argc <= 1) {
    std::cerr << "Required arguments not specified. Use '--help' for argument formats." << std::endl;
    std::cerr << help_description << std::endl;
    std::exit(EXIT_FAILURE);
  }

  auto parsed_opt = parseArgs(argc, argv);
  if (not parsed_opt) {
    fatalExit("Problem Parsing Command Line. Use '--help' for argument formats.");
  }

  auto const& parsed = *parsed_opt;

  if (parsed.helpRequested) {

    std::cerr << help_description << std::endl;
    std::exit(EXIT_SUCCESS);

  }

  // Work directory.
  if (parsed.workDirectory.empty()) {
    fatalExit("--workDirectory was not specified");
  }
  args_.workDirectory = parsed.workDirectory;
  std::cerr << "directory:" << args_.workDirectory << " was specified" << std::endl;

  bool valid_directory = Utility::directoryExists(getArgs().workDirectory);

  if (!valid_directory) {
    fatalExit("Specified work directory:" + getArgs().workDirectory + " does not exist.");
  }

  // Log file.
  if (parsed.logFile.empty()) {
    fatalExit("--logFile was not specified");
  }

  std::string log_file_name = parsed.logFile;
  // Join the log file and the directory
  std::string log_file_path = Utility::filePath(log_file_name, getArgs().workDirectory);
  // truncate the log file.
  std::fstream log_file(log_file_path, std::fstream::out | std::fstream::trunc);
  if (!log_file) {
    fatalExit("Cannot open log file (--logFile):" + log_file_path);
  }

  args_.logFile = log_file_path;

  // Options file.
  if (parsed.optionFile.empty()) {
    fatalExit("--optionFile was not specified");
  }
  args_.options_file = parsed.optionFile;

  return true;

}


/// Creates and returns the application logger.
std::unique_ptr<kel::ExecEnvLogger> kgl::GeneExecEnv::createLogger() {

  // Setup the Logger.
  return ExecEnv::createLogger(MODULE_NAME, getArgs().logFile, getArgs().max_error_count, getArgs().max_warn_count);

}