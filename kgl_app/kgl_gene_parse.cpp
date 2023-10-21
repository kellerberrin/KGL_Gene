//
// Created by kellerberrin on 4/05/18.
//

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Setup logger and read the XML program options.

#include "kgl_gene_app.h"
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include "kel_utility.h"

#include <iostream>
#include <fstream>

// Define namespace alias
namespace fs = boost::filesystem;
namespace po = boost::program_options;
namespace bt = boost;
namespace kgl = kellerberrin::genome;
namespace kel = kellerberrin;



// Parse the command line.
bool kgl::GeneExecEnv::parseCommandLine(int argc, char const ** argv)
{

  std::stringstream ss;
  ss << "Population Genome Comparison, module: "
     << MODULE_NAME
     << " version: "
     << VERSION << '\n'
     << "Usage: --workDirectory=<work_directory> --newLogFile=<new_log_file> --optionFile=<option_file.xml> (all arguments required)";
  const char* help_flag = "help";

  if (argc <= 1) {

    std::cerr << "Required arguments not specified. Use '--help' for argument formats." << std::endl;
    std::cerr << ss.str() << std::endl;
    std::exit(EXIT_FAILURE);

  }

  po::options_description runtime_options(ss.str());


  // Work directory
  const char* dir_desc =
      R"(The work directory where log files and data files are found.
  Use a Linux style directory specification with trailing forward slash '/'.
  Important - to run 'kgl_Gene' this directory must exist, it will not be created.)";
  const char* work_directory_flag = "workDirectory";

  // XML Options file
  const char* option_desc =
  R"(Specifies the options file (default "Runtime_Options.xml"). This file contains all runtime options.
     The file path is relative to the work directory)";
  const char* option_flag = "optionFile";

  // Logging file
  const char* log_desc =
  R"(Log file. Appends the log to any existing logs. The log file always resides in the work directory.)";
  const char* log_file_flag = "logFile";

  runtime_options.add_options ()
      (help_flag, ss.str().c_str())
      (work_directory_flag, po::value<std::string>(), dir_desc)
      (option_flag, po::value<std::string>(), option_desc)
      (log_file_flag, po::value<std::string>(), log_desc);

  po::variables_map variable_map;

  try {

    po::store(po::command_line_parser(argc, argv).options(runtime_options).run(), variable_map);
    po::notify(variable_map);

  } catch (po::error& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    std::cerr << "Problem Parsing Command Line. Use '--help' for argument formats." << std::endl;
    std::cerr << MODULE_NAME << " exits" << std::endl;
    std::exit(EXIT_FAILURE);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Move the following options to the options file.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  if (variable_map.count(help_flag)) {

    std::cerr << ss.str() << std::endl;
    std::exit(EXIT_SUCCESS);

  }

  if (variable_map.count(work_directory_flag)) {

    args_.workDirectory = variable_map[work_directory_flag].as<std::string>();
    std::cerr << "directory:" << args_.workDirectory << " was specified" << std::endl;

  } else {

    std::cerr << work_directory_flag << " was not specified" << std::endl;
    std::cerr << MODULE_NAME << " exits" << std::endl;
    std::exit(EXIT_FAILURE);

  }

  fs::path directory_path = fs::path(getArgs().workDirectory);

  boost::system::error_code error_code;
  bool valid_directory = fs::exists(directory_path, error_code);

  if (!valid_directory) {

    std::cerr << "Specified work directory:" << directory_path.string() << " does not exist." << std::endl;
    std::cerr << MODULE_NAME << " exits" << std::endl;
    std::exit(EXIT_FAILURE);

  }

  if (error_code.value() != boost::system::errc::success) {

    std::cerr << "Error verifying work directory:" << directory_path.string()
              << ", error was:" << error_code.message() << std::endl;
    std::cerr << MODULE_NAME << " exits" << std::endl;
    std::exit(EXIT_FAILURE);

  }

  // Get the log or newlog fields.
  if (variable_map.count(log_file_flag)) {

    std::string log_file_name;
    log_file_name = variable_map[log_file_flag].as<std::string>();
    // Join the log file and the directory
    fs::path log_file_path = directory_path / fs::path(log_file_name);
    // truncate the log file.
    std::fstream log_file(log_file_path.string(), std::fstream::out | std::fstream::trunc);
    if (!log_file) {

      std::cerr << "Cannot open log file (--logFile):" << log_file_path.string() << std::endl;
      std::cerr << MODULE_NAME << " exits" << std::endl;
      std::exit(EXIT_FAILURE);

    }

    args_.logFile = log_file_path.string();

  } else { // log file not specified - complain and exit.

    std::cerr << log_file_flag << " was not specified" << std::endl;
    std::cerr << MODULE_NAME << " exits" << std::endl;
    std::exit(EXIT_FAILURE);

  }

  // Read the options file.
  if (variable_map.count(option_flag)) {

    args_.options_file = variable_map[option_flag].as<std::string>();

  } else {

    std::cerr << option_flag << " was not specified" << std::endl;
    std::cerr << MODULE_NAME << " exits" << std::endl;
    std::exit(EXIT_FAILURE);

  }

  return true;

}


std::unique_ptr<kel::ExecEnvLogger> kgl::GeneExecEnv::createLogger() {

  // Setup the Logger.
  return ExecEnv::createLogger(MODULE_NAME, getArgs().logFile, getArgs().max_error_count, getArgs().max_warn_count);

}
