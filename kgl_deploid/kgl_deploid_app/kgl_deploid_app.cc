//
// Created by kellerberrin on 4/05/18.
//

#include "kgl_deploid_app.h"
#include "kgl_deploid_main.h"
#include "kgl_utility.h"

#include <iostream>
#include <cctype>
#include <seqan/arg_parse.h>
#define BOOST_FILESYSTEM_NO_DEPRECATED // Recommended by boost filesystem documentation.
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>


namespace kgl = kellerberrin::genome;
namespace fs = boost::filesystem;
namespace bt = boost;


// Static private member declarations.
kgl::DeploidArgs kgl::DeploidExecEnv::args_;

// Public static member functions.
const kgl::DeploidArgs& kgl::DeploidExecEnv::getArgs() { return args_; }

// Constants for the executable.
constexpr const char* kgl::DeploidExecEnv::MODULE_NAME;
constexpr const char* kgl::DeploidExecEnv::VERSION;


void kgl::DeploidExecEnv::executeApp() {

  const DeploidArgs& args = getArgs();

  kgl::deploidMain(args.argc, args.argv);

}


// Utility function to pre-pend work directory path and check file existance
void getFilePath(const std::string& option_text,
                 seqan::ArgumentParser& parser,
                 fs::path& directory_path,
                 std::string& file_name) {

  if (seqan::isSet(parser, option_text)) seqan::getOptionValue(file_name, parser, option_text);

  fs::path file_path = directory_path / fs::path(file_name);

  boost::system::error_code error_code;
  bool file_exists = fs::exists(file_path, error_code);

  if (!file_exists) {

    kgl::ExecEnv::log().critical("File: {}, type: {} does not exist.", file_path.string(), option_text);

  }

  if (error_code.value() != boost::system::errc::success) {

    kgl::ExecEnv::log().critical("File: {}, type: {}, error: {}."
        , file_path.string(), option_text, error_code.message());

  }

  file_name = file_path.string();

}


bool kgl::DeploidExecEnv::parseCommandLine(int argc, char const ** argv) {

  args_.argc = argc;
  args_.argv = argv;

  // Save the command line.
  ExecEnv::getCommandLine(argc, argv);

  // Setup ArgumentParser.
  seqan::ArgumentParser parser(MODULE_NAME);
  // Set short description, version, and date.
  setShortDescription(parser, "kgl_deploid COI haplotype generation");
  setVersion(parser, VERSION);
  setDate(parser, "May 2018");

  const char* program_desc =
      R"(kgl_deploid Complexity Of Infection (COI) haplotype generation)";

  addDescription(parser, program_desc);

  const char* dir_desc =
      R"(The work directory where log files and data files are found.
  Use a Linux style directory specification with trailing forward
  slash '/' (default './').
  Important - to run 'kgl_dEploid' this directory must exist, it will not be created.)";

  const char* workDirectoryFlag_ = "workDirectory";
  const char* workDirectoryShortFlag_ = "d";

  addOption(parser, seqan::ArgParseOption(workDirectoryShortFlag_, workDirectoryFlag_, dir_desc, seqan::ArgParseArgument::OUTPUT_DIRECTORY, "DIRECTORY"));

  const char* log_desc =
      R"(Log file. Appends the log to any existing logs (default "kgl_snp.log").'
  'The log file always resides in the work directory.)";

  const char* logFileFlag_ = "logFile";
  const char* logFileShortFlag_ = "l";

  addOption(parser, seqan::ArgParseOption(logFileShortFlag_, logFileFlag_, log_desc, seqan::ArgParseArgument::OUTPUT_FILE, "LOG_FILE"));

  const char* newlog_desc =
      R"(Flush an existing log file (file name argument optional, (default "kgl_snp.log").
  The log file always resides in the work directory.)";

  const char* newLogFileFlag_ = "newLogFile";
  const char* newLogFileShortFlag_ = "n";

  addOption(parser, seqan::ArgParseOption(newLogFileShortFlag_, newLogFileFlag_, newlog_desc, seqan::ArgParseArgument::OUTPUT_FILE, "TRUNC_LOG_FILE"));

  const char* vcf_file_desc =
  R"(The VCF file path. This is relative to the work directory. and the workDirectory is prepended to to this path e.g. 'workDirectory/vcf_file_path')";

  const char* vcfFileFlag_ = "vcfFile";
  const char* vcfFileShortFlag_ = "vcf";

  addOption(parser, seqan::ArgParseOption(vcfFileShortFlag_, vcfFileFlag_, vcf_file_desc, seqan::ArgParseArgument::INPUT_FILE, "VCF_FILE"));

  const char* plaf_file_desc =
  R"(The Population Allele Frequency file path (plaf). This is relative to the work directory. and the workDirectory
  is prepended to to this path e.g. 'workDirectory/plaf_file_path')";

  const char* plafFileFlag_ = "plafFile";
  const char* plafFileShortFlag_ = "plaf";

  addOption(parser, seqan::ArgParseOption(plafFileShortFlag_, plafFileFlag_, plaf_file_desc, seqan::ArgParseArgument::INPUT_FILE, "INPUT_FILE"));



  // Parse command line.
  seqan::ArgumentParser::ParseResult parse_result = seqan::parse(parser, argc, argv);

  if (parse_result == seqan::ArgumentParser::PARSE_HELP) {

    std::exit(EXIT_SUCCESS);  // Help, -h or --help flags used, so just exit.

  } else if (parse_result != seqan::ArgumentParser::PARSE_OK) {  // Problem parsing the command line.

    std::cerr << "Error - Problem Parsing Command Line. Use '-h' or '--help' for argument formats." << std::endl;
    std::cerr << MODULE_NAME << " exits" << std::endl;
    std::exit(EXIT_FAILURE);

  }

  // Get the work directory and check that it exists
  if (seqan::isSet(parser, workDirectoryFlag_)) seqan::getOptionValue(args_.workDirectory, parser, workDirectoryFlag_);

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
  if (seqan::isSet(parser, newLogFileFlag_)) {

    std::string log_file_name;
    seqan::getOptionValue(log_file_name, parser, newLogFileFlag_);
    // Join the log file and the directory
    fs::path log_file_path = directory_path / fs::path(log_file_name);
    // truncate the log file.
    std::fstream log_file(log_file_path.string(), std::fstream::out | std::fstream::trunc);
    if(!log_file)
    {
      std::cerr << "Cannot open truncated log file (--newLogFile):" << log_file_path.string() << std::endl;
      std::cerr << MODULE_NAME << " exits" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    args_.logFile = log_file_path.string();

  } else if (seqan::isSet(parser, logFileFlag_)) {

    std::string log_file_name;
    seqan::getOptionValue(log_file_name, parser, logFileFlag_);
    // Join the log file and the directory
    fs::path log_file_path = directory_path / fs::path(log_file_name);
    args_.logFile = log_file_path.string();

  } else { // log file not specified - join default log file to the work directory.

    fs::path log_file_path = directory_path / fs::path(args_.logFile);
    args_.logFile = log_file_path.string();

  }
  // Setup the Logger.
  ExecEnv::createLogger(MODULE_NAME, getArgs().logFile, getArgs().max_error_count, getArgs().max_warn_count);

  return true;

}