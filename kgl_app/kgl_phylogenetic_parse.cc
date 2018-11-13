//
// Created by kellerberrin on 4/05/18.
//



#include <iostream>
#include <cctype>
#include <seqan/arg_parse.h>
#include "kgl_phylogenetic_app.h"
#include "kgl_phylogenetic_app_analysis.h"
#define BOOST_FILESYSTEM_NO_DEPRECATED // Recommended by boost filesystem documentation.
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include "kgl_utility.h"


// Define namespace alias
namespace fs = boost::filesystem;
namespace bt = boost;
namespace kgl = kellerberrin::genome;


// Static private member declarations.
kgl::Phylogenetic kgl::PhylogeneticExecEnv::args_;
// Static storage for the runtime options.
kgl::RuntimeProperties kgl::PhylogeneticExecEnv::runtime_options_;

// Public static member functions.
const kgl::Phylogenetic& kgl::PhylogeneticExecEnv::getArgs() { return args_; }
const kgl::RuntimeProperties& kgl::PhylogeneticExecEnv::getRuntimeOptions() { return runtime_options_; }

// Constants for the executable.
constexpr const char* kgl::PhylogeneticExecEnv::MODULE_NAME;
constexpr const char* kgl::PhylogeneticExecEnv::VERSION;



// Parse the command line.
bool kgl::PhylogeneticExecEnv::parseCommandLine(int argc, char const ** argv)
{
  // Save the command line.
  ExecEnv::getCommandLine(argc, argv);

  // Setup ArgumentParser.
  seqan::ArgumentParser parser(MODULE_NAME);
  // Set short description, version, and date.
  setShortDescription(parser, "Population Genome Comparison");
  setVersion(parser, VERSION);
  setDate(parser, "September 2018");

  // Define usage line and long description.
  addUsageLine(parser, "--workDirectory <work.directory>  --newLogFile <new_log_file> --optionFile <optionFile>");

  const char* program_desc =
      R"("kgl_phylogenetic" is a fast C++ program to find genetic differences (SNPs/Indels) in a population of organisms.
  The entire genome of many organisms can be compared and analysed simultaneously (with sufficient memory).
  The --fileList flag" specifies a list of SAM/BAM/VCF files and associated organism names to be processed.
  This program also takes the genome FASTA file and the corresponding genetic feature  model in GFF3 (only) format
  and builds a memory database of the genetic structure of the target organism to facilitate analysis.
  Source FASTQ files should be filtered for quality and aligned to the FASTA reference genome with a tool
  such as BowTie2 or bwa, and output in SAM/BAM format. The program will read VCF files without further processing)";

  addDescription(parser, program_desc);

  // Define Options -- Section Modification Options
  addSection(parser, "Required Program Options");

  const char* dir_desc =
      R"(The work directory where log files and data files are found.
  Use a Linux style directory specification with trailing forward
  slash '/' (default './Work/').
  Important - to run 'kgl_snp' this directory must exist, it will not be created.)";

  const char* workDirectoryFlag_ = "workDirectory";
  const char* workDirectoryShortFlag_ = "d";

  addOption(parser, seqan::ArgParseOption(workDirectoryShortFlag_, workDirectoryFlag_, dir_desc, seqan::ArgParseArgument::OUTPUT_DIRECTORY, "DIRECTORY"));

  const char* option_desc =
  R"(Specifies the options file. This file contains all runtime options. The file path is relative to the work directory)";

  const char* optionFlag_ = "optionFile";
  const char* optionShortFlag_ = "e";

  addOption(parser, seqan::ArgParseOption(optionShortFlag_, optionFlag_, option_desc, seqan::ArgParseArgument::INPUT_FILE, "OPTION_FILE"));

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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Move the following options to the options file.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  const char* analysis_desc =
      R"(Define which genome analytic to perform. Defaults to '*' for all analytics.)";

  const char* analysisFlag_ = "analysisType";
  const char* analysisShortFlag_ = "z";

  addOption(parser, seqan::ArgParseOption(analysisShortFlag_, analysisFlag_, analysis_desc, seqan::ArgParseArgument::STRING, "ANALYSIS_TYPE"));

  const char* outCSV_desc =
      R"(Input an comma delimiter CSV file for further, analysis specific, processing, (default "/...WorkDirectory.../kgl_aux.csv").
  This file is always relative to the work directory.)";

  const char* CSVFileFlag_ = "auxCSVFile";
  const char* CSVFileShortFlag_ = "o";

  addOption(parser, seqan::ArgParseOption(CSVFileShortFlag_, CSVFileFlag_, outCSV_desc, seqan::ArgParseArgument::OUTPUT_FILE, "AUX_CSV_FILE"));

  const char* verbose_desc =
      R"(Flag. Enables verbose logged output to screen and log file.)";

  const char* verboseFlag_ = "verbose";
  const char* verboseShortFlag_ = "v";

  addOption(parser, seqan::ArgParseOption(verboseShortFlag_, verboseFlag_, verbose_desc, seqan::ArgParseArgument::BOOL, "FLAG"));

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

  // Read the options file.
  if (seqan::isSet(parser, optionFlag_)) {

    std::string options_file;
    seqan::getOptionValue(options_file, parser, optionFlag_);
    std::string options_file_path = kgl::Utility::filePath(options_file, getArgs().workDirectory);
    if (not runtime_options_.readProperties(options_file_path)) {

      ExecEnv::log().critical("parseCommandLine; could not read specified runtime properties file: {}", options_file_path);

    }

  }


  if (seqan::isSet(parser, analysisFlag_)) seqan::getOptionValue(args_.analysisType, parser, analysisFlag_);

  args_.analysisType = kgl::Utility::toupper(args_.analysisType);

  if (not PhylogeneticAnalysis::checkAnalysisType(args_.analysisType)) {

    ExecEnv::log().critical("Invalid Analysis Requested: {}; program terminates.", args_.analysisType);

  }

  if (seqan::isSet(parser, verboseFlag_)) seqan::getOptionValue(args_.verbose, parser, verboseFlag_);

  ExecEnv::log().SetVerbose(getArgs().verbose);  // Set the logger verbose level.

  return true;

}

