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


// Simple command-line parser.
// Supports --flag=value and --flag value syntax for three required string flags
// plus --help/-h as a boolean flag.
struct ParsedArgs {

  std::string reference_directory;
  std::string log_file;
  std::string option_file;
  std::string option_file_out;
  std::string parsed_option_out;
  bool help_requested{false};

};


/// Parses command-line arguments into a ParsedArgs structure, supporting --flag=value and --flag value syntax.
[[nodiscard]] static std::optional<ParsedArgs> parseArgs(int argc, char const** argv) {

  ParsedArgs args;

  // Map flag names to the ParsedArgs members they populate.
  static const std::unordered_map<std::string, std::string*> string_flags = {
    {"--referenceDirectory", &args.reference_directory},
    {"--logFile",       &args.log_file},
    {"--optionFile",    &args.option_file},
    {"--optionFileOut",    &args.option_file_out},
    {"--parsedTree",    &args.parsed_option_out},
  };

  for (int i = 1; i < argc; ++i) {

    std::string arg = argv[i];

    if (arg == "--help" or arg == "-h") {

      args.help_requested = true;
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

kgl::CmdLineArgs validateArguments(const ParsedArgs& parsed_args) {

  // Reference directory non empty
  if (parsed_args.reference_directory.empty()) {

    fatalExit("--workDirectory was not specified");

  }

  // Reference directory exists
  bool valid_directory = kel::Utility::directoryExists(parsed_args.reference_directory);
  if (!valid_directory) {

    fatalExit("Specified work directory:" + parsed_args.reference_directory + " does not exist.");

  }

  // Log file non empty
  if (parsed_args.log_file.empty()) {

    fatalExit("--logFile was not specified");

  }

  // Options file not empty
  if (parsed_args.option_file.empty()) {

    fatalExit("--optionFile was not specified");

  }

  if (parsed_args.option_file == parsed_args.log_file) {

    fatalExit("--optionFile cannot have the same file name as --logFile.");

  }


  // Check the optional arguments.
  if (not parsed_args.option_file_out.empty()) {

    if (parsed_args.option_file_out == parsed_args.option_file) {

      fatalExit("--optionFileOut cannot have the same file name as --optionFile.");

    }

    if (parsed_args.option_file_out == parsed_args.log_file) {

      fatalExit("--optionFileOut cannot have the same file name as --logFile.");

    }

    if (not parsed_args.option_file_out.empty() and not parsed_args.parsed_option_out.empty()) {

      if (parsed_args.option_file_out == parsed_args.parsed_option_out) {

        fatalExit("--optionFileOut cannot have the same file name as --parsedTree");

      }

    }

  }

  if (not parsed_args.parsed_option_out.empty()) {

    if (parsed_args.parsed_option_out == parsed_args.option_file) {

      fatalExit("--parsedTree cannot have the same file name as --optionFile.");

    }

    if (parsed_args.option_file_out == parsed_args.log_file) {

      fatalExit("--parsedTree cannot have the same file name as --logFile.");

    }

  }

  kgl::CmdLineArgs cmd_args;

  cmd_args.reference_directory = parsed_args.reference_directory;
  cmd_args.log_file = kel::Utility::filePath(parsed_args.log_file, cmd_args.reference_directory);
  cmd_args.options_file = kel::Utility::filePath(parsed_args.option_file, cmd_args.reference_directory);
  cmd_args.option_file_out = kel::Utility::filePath(parsed_args.option_file_out, cmd_args.reference_directory);
  cmd_args.parsed_option_out = kel::Utility::filePath(parsed_args.parsed_option_out, cmd_args.reference_directory);

  return cmd_args;

}

/// Parses the command line arguments and initializes the runtime environment.
bool kgl::GeneExecEnv::parseCommandLine(int argc, char const ** argv)
{

  std::stringstream ss;
  ss << "Population Genome Comparison, module: "
     << MODULE_NAME
     << " version: "
     << VERSION << '\n'
     << "Required Arguments: --workDirectory=<work_directory> --logFile=<log_file.log> --optionFile=<option_file.xml>" << '\n'
     << "Optional Arguments: --optionFileOut=<option_file_out.xml> --parsedTree=<parsed_tree.txt>" << '\n'
     << "The optional arguments are used for debugging the option file and parsed option tree." << '\n'
     << "To prevent accidental file overwrite it is a requirement that:" << '\n'
     <<  "<option_file.xml> != <option_file_out.xml>, <option_file.xml> != <parsed_tree.txt>, <option_file_out.xml> != <parsed_tree.txt>";

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

  if (parsed.help_requested) {

    std::cerr << help_description << std::endl;
    std::exit(EXIT_SUCCESS);

  }

  args_ = validateArguments(parsed);

  // truncate the log file.
  std::fstream log_file(args_.log_file, std::fstream::out | std::fstream::trunc);
  if (!log_file) {

    fatalExit("Cannot open log file (--logFile):" + args_.log_file);

  }

  return true;

}

/// Creates and returns the application logger.
std::unique_ptr<kel::ExecEnvLogger> kgl::GeneExecEnv::createLogger() {

  // Setup the Logger.
  return ExecEnv::createLogger(MODULE_NAME, getArgs().log_file, getArgs().max_error_count, getArgs().max_warn_count);

}


