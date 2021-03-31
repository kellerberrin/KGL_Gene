//
// Created by kellerberrin on 4/05/18.
//

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Currently uses the Seqan library to parse arguments. Switch to boost when convenient.
// Seqan dependencies should be limited to sequence analysis only.

#include <iostream>
#include <seqan/arg_parse.h>
#include "kgl_gene_app.h"
#include <boost/filesystem.hpp>
#include "kel_utility.h"


// Define namespace alias
namespace fs = boost::filesystem;
namespace bt = boost;
namespace kgl = kellerberrin::genome;



// Public static member functions.

// Parse the command line.
bool kgl::GeneExecEnv::parseCommandLine(int argc, char const ** argv)
{

  // Setup ArgumentParser.
  seqan::ArgumentParser parser(MODULE_NAME);
  // Set short description, version, and date.
  setShortDescription(parser, "Population Genome Comparison");
  setVersion(parser, VERSION);
  setDate(parser, "March 2021");

  // Define usage line and long description.
  addUsageLine(parser, "--workDirectory <work.directory>  --newLogFile <new_log_file> --optionFile <optionFile>");

  const char* program_desc =
      R"("kgl_Gene" analyzes genetic differences (SNPs/Indels) in a population of organisms.
  The entire genome of many organisms can be compared and analysed simultaneously (with sufficient memory).
  The options xml file specifies a list of VCF files and associated organism attributes to be processed.
  This program also takes the genome FASTA file and the corresponding genetic feature model in GFF3 (only) format
  and builds a memory database of the genetic structure of the target organism to facilitate analysis.)";

  addDescription(parser, program_desc);

  // Define Options -- Section Modification Options
  addSection(parser, "Command Line Program Options");

  const char* dir_desc =
      R"(The work directory where log files and data files are found.
  Use a Linux style directory specification with trailing forward
  slash '/' (default './Work/').
  Important - to run 'kgl_snp' this directory must exist, it will not be created.)";

  const char* workDirectoryFlag_ = "workDirectory";
  const char* workDirectoryShortFlag_ = "d";

  addOption(parser, seqan::ArgParseOption(workDirectoryShortFlag_, workDirectoryFlag_, dir_desc, seqan::ArgParseArgument::OUTPUT_DIRECTORY, "DIRECTORY"));

  const char* option_desc =
  R"(Specifies the options file (default "Runtime_Options.xml"). This file contains all runtime options. The file path is relative to the work directory)";

  const char* optionFlag_ = "optionFile";
  const char* optionShortFlag_ = "o";

  addOption(parser, seqan::ArgParseOption(optionShortFlag_, optionFlag_, option_desc, seqan::ArgParseArgument::INPUT_FILE, "OPTION_FILE"));

  const char* log_desc =
  R"(Log file. Appends the log to any existing logs (default "kgl_snp.log").'
  'The log file always resides in the work directory.)";

  const char* logFileFlag_ = "logFile";
  const char* logFileShortFlag_ = "l";

  addOption(parser, seqan::ArgParseOption(logFileShortFlag_, logFileFlag_, log_desc, seqan::ArgParseArgument::OUTPUT_FILE, "LOG_FILE"));


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Move the following options to the options file.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
  if (seqan::isSet(parser, logFileFlag_)) {

    std::string log_file_name;
    seqan::getOptionValue(log_file_name, parser, logFileFlag_);
    // Join the log file and the directory
    fs::path log_file_path = directory_path / fs::path(log_file_name);
    // truncate the log file.
    std::fstream log_file(log_file_path.string(), std::fstream::out | std::fstream::trunc);
    if(!log_file)
    {
      std::cerr << "Cannot open log file (--logFile):" << log_file_path.string() << std::endl;
      std::cerr << MODULE_NAME << " exits" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    args_.logFile = log_file_path.string();

  } else { // log file not specified - join default log file to the work directory.

    fs::path log_file_path = directory_path / fs::path(args_.logFile);
    args_.logFile = log_file_path.string();

  }
  // Setup the Logger.
  ExecEnv::createLogger(MODULE_NAME, getArgs().logFile, getArgs().max_error_count, getArgs().max_warn_count);

  // Read the options file.
  if (seqan::isSet(parser, optionFlag_)) {

    std::string options_file;
    seqan::getOptionValue(options_file, parser, optionFlag_);
    args_.options_file = options_file;
  }

  runtime_options_.setWorkDirectory(args_.workDirectory);
  if (not runtime_options_.readProperties(args_.options_file)) {

    std::string options_file_path = Utility::filePath(args_.options_file, args_.workDirectory);
    ExecEnv::log().critical("parseCommandLine; could not read specified runtime properties file: {}", options_file_path);

  }


  return true;

}

