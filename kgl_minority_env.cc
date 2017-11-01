///
// Created by kellerberrin on 17/10/17.
//


#include <iostream>
#include <seqan/arg_parse.h>
#include "kgl_minority_env.h"
#define BOOST_FILESYSTEM_NO_DEPRECATED // Recommended by boost filesystem documentation.
#include <boost/filesystem.hpp>


// Define namespace alias
namespace fs = boost::filesystem;
namespace kgl = kellerberrin::genome;


// Static private member declarations.
kgl::MinorityArgs kgl::MinorityExecEnv::args_;

// Public static member functions.
const kgl::MinorityArgs& kgl::MinorityExecEnv::args() { return args_; }

// Constants for the executable.
constexpr const char* kgl::MinorityExecEnv::MODULE_NAME;
constexpr const char* kgl::MinorityExecEnv::VERSION;


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

    kgl::ExecEnv::log().critical("File: {}, type: {} does not exist."
    , file_path.string(), option_text);

  }

  if (error_code.value() != boost::system::errc::success) {

    kgl::ExecEnv::log().critical("File: {}, type: {}, error: {}."
    , file_path.string(), option_text, error_code.message());

  }

  file_name = file_path.string();

}


bool readFileList(const std::string& file_list_name,
                  fs::path& directory_path,
                  std::vector<std::string>& file_list) {

  bool result = true;
  std::ifstream list_file;

  file_list.clear();
  // Open input file.
  list_file.open(file_list_name);

  if (not list_file.good()) {

    return false;

  }

  try {

    while (true) {

      std::string sam_file_name;

      if (std::getline(list_file, sam_file_name).eof()) break;

      fs::path file_path = directory_path / fs::path(sam_file_name);

      boost::system::error_code error_code;
      bool file_exists = fs::exists(file_path, error_code);

      if (!file_exists) {

        kgl::ExecEnv::log().critical("File: {} does not exist.", file_path.string());

      }

      if (error_code.value() != boost::system::errc::success) {

        kgl::ExecEnv::log().critical("File: {}, error: {}.", file_path.string(), error_code.message());

      }

      file_list.push_back(file_path.string());

    }

    list_file.close();

  }
  catch (std::exception const &e) {

    result = false;

  }

  return result;

}

// Parse the command line.
bool kgl::MinorityExecEnv::parseCommandLine(int argc, char const ** argv)
{
  // Setup ArgumentParser.
  seqan::ArgumentParser parser(MODULE_NAME);
  // Set short description, version, and date.
  setShortDescription(parser, "Haploid Organism Genome Comparison");
  setVersion(parser, VERSION);
  setDate(parser, "October 2017");

  // Define usage line and long description.
  addUsageLine(parser, "-d <work.directory>  -f <ref.fasta> -g <ref.gff> -m <mutant.sam>");

  const char* program_desc =
  R"("kgl_minority" is a very fast C++ program to find genetic differences (SNPs) in (optional)
  parent and mutant pairs by comparing genomic sequencing reads aligned to the same reference genome.
  If the optional --parent flag is specified then the mutant genome is compared to the parent genome and
  the reference genome. The parent genome is also compared to the reference genome. If the --parent flag"
  is not specified then the --mutant is compared to the reference FASTA genome only.
  This program takes the genome FASTA file, the corresponding gene model in GFF3 (only) format, and the
  (optional) parent and mutant SAM files, respectively. Source FASTQ files should be filtered for quality
  and aligned to the FASTA reference genome with a tool such as BowTie2 or bwa, and output in SAM/BAM format.)";

  addDescription(parser, program_desc);


  // Define Options -- Section Modification Options
  addSection(parser, "Modification Options");

  const char* dir_desc =
  R"(The work directory where log files and data files are found.
  Use a Linux style directory specification with trailing forward
  slash '/' (default './Work/').
  Important - to run 'kgl_snp' this directory must exist, it will not be created.)";

  addOption(parser, seqan::ArgParseOption("d", "workDirectory", dir_desc,
                                          seqan::ArgParseArgument::OUTPUT_DIRECTORY, "DIRECTORY"));

  const char* gff_desc =
  R"(The gff3 (not GFF2 or GTF) file that contains the genes, exons, etc for the
  chromosome(s)/contiguous region(s) of interest)";

  addOption(parser, seqan::ArgParseOption("g", "gffFile", gff_desc,
                                          seqan::ArgParseArgument::INPUT_FILE, "GFF_FILE"));

  const char* fasta_desc =
  R"(The fasta file that contains the reference nucleotide sequence for the chromsome(s) of interest)";

  addOption(parser, seqan::ArgParseOption("f", "fastaFile", fasta_desc,
                                          seqan::ArgParseArgument::INPUT_FILE, "FASTA_FILE"));

  const char* parent_desc =
  R"(The SAM file generated by sequencing the parent organism. The SAM file reads should be
   aligned with respect to the reference genome specified by the --fasta flag.)";

  addOption(parser, seqan::ArgParseOption("p", "parentFile", parent_desc,
                                          seqan::ArgParseArgument::INPUT_FILE, "PARENT_READS"));

  const char* mutant_desc =
  R"(The SAM file generated by sequencing the MUTANT organism. The SAM file reads should be
   aligned with respect to the reference genome specified by the --fasta flag.)";

  addOption(parser, seqan::ArgParseOption("m", "mutantFile", mutant_desc,
                                          seqan::ArgParseArgument::INPUT_FILE, "MUTANT_READS"));

  const char* file_list_desc =
  R"(A file with a list of all SAM files to be read and analyzed. The SAM file reads should be
   aligned with respect to the reference genome specified by the --fasta flag.)";

  addOption(parser, seqan::ArgParseOption("z", "fileList", file_list_desc,
                                          seqan::ArgParseArgument::INPUT_FILE, "FILE_LIST_READS"));

  const char* contig_desc =
  R"(Define which contiguous DNA region (chromosome/mitochondria) to process.
  Defaults to '*' for all contiguous regions.)";

  addOption(parser, seqan::ArgParseOption("c", "contig", contig_desc,
                                          seqan::ArgParseArgument::STRING, "CONTIG_ID"));

  const char* log_desc =
  R"(Log file. Appends the log to any existing logs (default "kgl_snp.log").'
  'The log file always resides in the work directory.)";

  addOption(parser, seqan::ArgParseOption("l", "logFile", log_desc,
                                          seqan::ArgParseArgument::OUTPUT_FILE, "LOG_FILE"));

  const char* newlog_desc =
  R"(Flush an existing log file (file name argument optional, (default "kgl_snp.log").
  The log file always resides in the work directory.)";

  addOption(parser, seqan::ArgParseOption("n", "newLogFile", newlog_desc,
                                          seqan::ArgParseArgument::OUTPUT_FILE, "TRUNC_LOG_FILE"));

  const char* mutant_min_count_desc =
  R"(The minimum SAM/BAM read coverage for the a single nucleotide analyzed in the Mutant genome.)";

  addOption(parser, seqan::ArgParseOption("mmc", "mutantMinimumCount", mutant_min_count_desc,
                                          seqan::ArgParseArgument::INTEGER, "INT"));

  const char* mutant_min_prop_desc =
  R"(The minimum proportion of a single nucleotide analyzed in the Mutant genome that is
  at variance from the reference (fasta) genome)";

  addOption(parser, seqan::ArgParseOption("mmp", "mutantMinimumProportion", mutant_min_prop_desc,
                                          seqan::ArgParseArgument::DOUBLE, "FLOAT"));

  const char* parent_min_count_desc =
  R"(The minimum SAM/BAM count of a single nucleotide analyzed in the Parent (wild-type)
  genome that is at variance from the reference (fasta) genome)";

  addOption(parser, seqan::ArgParseOption("pmc", "parentMinimumCount", parent_min_count_desc ,
                                          seqan::ArgParseArgument::INTEGER, "INT"));

  const char* parent_min_prop_desc =
  R"(The minimum proportion of a single nucleotide analyzed in the Parent (wild-type)
  genome that is at variance from the reference (fasta) genome)";

  addOption(parser, seqan::ArgParseOption("pmp", "parentMinimumProportion", parent_min_prop_desc,
                                          seqan::ArgParseArgument::DOUBLE, "FLOAT"));

  const char* thread_count_desc =
  R"(The number of CPU processes/threads assigned to processing genome data.
  Defaults to the number of CPU processors available (-1).)";

  addOption(parser, seqan::ArgParseOption("t", "threadCount", thread_count_desc,
                                          seqan::ArgParseArgument::INTEGER, "INT"));

  const char* queue_size_desc =
  R"(The maximum number of SAM/BAM records held in the inter-process record queue (default 1000000).)";

  addOption(parser, seqan::ArgParseOption("q", "queueSize", queue_size_desc,
                                          seqan::ArgParseArgument::INTEGER, "INT"));

  const char* translation_table_desc =
  R"(The amino acid translation table used in gene verification and amino acid generation from CDS regions.
  The table code is an integer that matches the translation tables defines on the NCBI website at
  https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)";

  addOption(parser, seqan::ArgParseOption("a", "translationTable", translation_table_desc,
                                          seqan::ArgParseArgument::INTEGER, "INT"));

  const char* verbose_desc =
  R"(Flag. Enables verbose logged output to screen and log file.)";

  addOption(parser, seqan::ArgParseOption("v", "verbose", verbose_desc,
                                          seqan::ArgParseArgument::BOOL, "FLAG"));

  const char* read_quality_desc =
  R"(The nucleotide read quality as -10 log10 Pr {ReadError} e.g. 30 is a 1/1000 chance
  of an nucleotide read error. Defaults to 30. Set this value to 0 to disable quality checking.)";

  addOption(parser, seqan::ArgParseOption("rq", "readQuality", read_quality_desc,
                                          seqan::ArgParseArgument::INTEGER, "INT"));

  const char* lock_granularity_desc =
  R"(The number of nucleotide positions per inter-process write lock (less is faster).)";

  addOption(parser, seqan::ArgParseOption("lg", "lockGranularity", lock_granularity_desc,
                                          seqan::ArgParseArgument::INTEGER, "INT"));


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
  if (seqan::isSet(parser, "workDirectory")) seqan::getOptionValue(args_.workDirectory, parser, "workDirectory");

  fs::path directory_path = fs::path(args().workDirectory);

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
  if (seqan::isSet(parser, "newLogFile")) {

    std::string log_file_name;
    seqan::getOptionValue(log_file_name, parser, "newLogFile");
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

  } else if (seqan::isSet(parser, "logFile")) {

    std::string log_file_name;
    seqan::getOptionValue(log_file_name, parser, "logFile");
    // Join the log file and the directory
    fs::path log_file_path = directory_path / fs::path(log_file_name);
    args_.logFile = log_file_path.string();

  } else { // log file not specified - join default log file to the work directory.

    fs::path log_file_path = directory_path / fs::path(args_.logFile);
    args_.logFile = log_file_path.string();

  }
  // Setup the Logger.
  createLogger(MODULE_NAME, args().logFile, args().max_error_count, args().max_warn_count);

  std::string file_list;
  // Setup the files.
  getFilePath("fastaFile" ,parser , directory_path, args_.fastaFile);
  getFilePath("gffFile" ,parser , directory_path, args_.gffFile);
  if (seqan::isSet(parser, "mutantFile")) getFilePath("mutantFile" ,parser , directory_path, args_.mutantFile);
  if (seqan::isSet(parser, "parentFile")) getFilePath("parentFile" ,parser , directory_path, args_.parentFile);
  if (seqan::isSet(parser, "fileList")) getFilePath("fileList" ,parser , directory_path, file_list);

  if (args_.mutantFile.length() == 0 and file_list.length() == 0) {

    ExecEnv::log().critical("No SAM files specified. Specify a SAM file using '--mutantFile' or a SAM list '--fileList'.");

  }

  if (args_.mutantFile.length() > 0 and file_list.length() > 0) {

    ExecEnv::log().warn("Both '--mutantFile' and '--fileList' specified; processing the '--fileList' of SAM files.");

  }

  if (file_list.length() > 0) {

    if (not readFileList(file_list, directory_path, args_.fileList)) {

      ExecEnv::log().critical("Unable to read specified file list: {}", file_list);

    }

  }

  // Setup Other Parameters.
  if (seqan::isSet(parser, "contig")) seqan::getOptionValue(args_.contig, parser, "contig");
  if (seqan::isSet(parser, "mutantMinimumCount")) seqan::getOptionValue(args_.mutantMinCount, parser,
                                                                        "mutantMinimumCount");
  if (seqan::isSet(parser, "mutantMinimumProportion")) seqan::getOptionValue(args_.mutantMinProportion, parser,
                                                                             "mutantMinimumProportion");
  if (seqan::isSet(parser, "parentMinimumCount")) seqan::getOptionValue(args_.mutantMinCount, parser,
                                                                        "parentMinimumCount");
  if (seqan::isSet(parser, "mutantMinimumProportion")) seqan::getOptionValue(args_.mutantMinProportion, parser,
                                                                             "parentMinimumProportion");
  if (seqan::isSet(parser, "queueSize")) seqan::getOptionValue(args_.queueSize, parser,"queueSize");
  if (seqan::isSet(parser, "threadCount")) seqan::getOptionValue(args_.queueSize, parser,"threadCount");
  if (seqan::isSet(parser, "lockGranularity")) seqan::getOptionValue(args_.queueSize, parser,"lockGranularity");
  if (seqan::isSet(parser, "translationTable")) seqan::getOptionValue(args_.aminoTranslationTable, parser,
                                                                      "translationTable");
  if (seqan::isSet(parser, "verbose")) seqan::getOptionValue(args_.verbose, parser, "verbose");

  if (seqan::isSet(parser, "readQuality")) {

    int readQuality;
    seqan::getOptionValue(readQuality, parser, "readQuality");
    args_.readQuality = static_cast<unsigned char>(readQuality);

  }

  return true;

}

