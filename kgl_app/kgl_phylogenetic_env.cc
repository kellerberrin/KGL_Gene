//
// Created by kellerberrin on 10/11/17.
//


#include <iostream>
#include <cctype>
#include <seqan/arg_parse.h>
#include "kgl_phylogenetic_env.h"
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

// Public static member functions.
const kgl::Phylogenetic& kgl::PhylogeneticExecEnv::args() { return args_; }

// Constants for the executable.
constexpr const char* kgl::PhylogeneticExecEnv::MODULE_NAME;
constexpr const char* kgl::PhylogeneticExecEnv::VERSION;


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
                  std::vector<kgl::FileListInfo>& file_list) {

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

      std::string line;
      std::vector<std::string> tokens;
      std::string sam_file_name;
      std::string genome_name;
      int line_number = 0;

      if (std::getline(list_file, line).eof()) break;

      line_number++;

      if (line.empty()) {  // Skip empty lines.

        continue;

      } else {

        if (line[0] == '#') {  // Skip comments.

          continue;

        }

      }

      boost::char_separator<char> sep("\t ");
      using tokenizer_t = bt::tokenizer< boost::char_separator<char>>;
      tokenizer_t tok(line, sep); // looking for 1 or 2 tokens separated by whitespace.
      for(auto it = tok.begin(); it != tok.end(); ++it){

        tokens.push_back(*it);

      }

      if (tokens.empty()) {

        kgl::ExecEnv::log().critical("No SAM/BAM filename given on line number: {} in list file: {}.",
                                     line_number, file_list_name);

      } else if (tokens.size() == 1) {

        sam_file_name = tokens[0];
        genome_name = tokens[0];

      } else {

        sam_file_name = tokens[0];
        genome_name = tokens[1];

      }

      fs::path file_path = directory_path / fs::path(sam_file_name);

      boost::system::error_code error_code;
      bool file_exists = fs::exists(file_path, error_code);

      if (!file_exists) {

        kgl::ExecEnv::log().critical("File: {} does not exist.", file_path.string());

      }

      if (error_code.value() != boost::system::errc::success) {

        kgl::ExecEnv::log().critical("File: {}, error: {}.", file_path.string(), error_code.message());

      }

      kgl::FileListInfo sam_bam_info;
      sam_bam_info.file_name = file_path.string();
      sam_bam_info.genome_name = genome_name;
      file_list.push_back(sam_bam_info);

    }

    list_file.close();

  }
  catch (std::exception const &e) {

    result = false;

  }

  return result;

}

std::string getCommandLine(int argc, char const ** argv) {

  std::string command_line;

  for (int idx = 0; idx < argc; ++idx) {

    command_line += argv[idx];
    command_line += ' ';

  }

  return command_line;

}


// Parse the command line.
bool kgl::PhylogeneticExecEnv::parseCommandLine(int argc, char const ** argv)
{
  // Get the command line.
  kgl::ExecEnv::commandLine(getCommandLine(argc, argv));

  // Setup ArgumentParser.
  seqan::ArgumentParser parser(MODULE_NAME);
  // Set short description, version, and date.
  setShortDescription(parser, "Population Genome Comparison");
  setVersion(parser, VERSION);
  setDate(parser, "December 2017");

  // Define usage line and long description.
  addUsageLine(parser, "-d <work.directory>  -f <ref.fasta> -g <ref.gff> -m <fileList>");

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
  addSection(parser, "Modification Options");

  const char* dir_desc =
  R"(The work directory where log files and data files are found.
  Use a Linux style directory specification with trailing forward
  slash '/' (default './Work/').
  Important - to run 'kgl_snp' this directory must exist, it will not be created.)";

  const char* workDirectoryFlag_ = "workDirectory";
  const char* workDirectoryShortFlag_ = "d";

  addOption(parser, seqan::ArgParseOption(workDirectoryShortFlag_, workDirectoryFlag_, dir_desc, seqan::ArgParseArgument::OUTPUT_DIRECTORY, "DIRECTORY"));

  const char* gff_desc =
  R"(The gff3 (not GFF2 or GTF) file that contains the genes, exons, etc for the
  chromosome(s)/contiguous region(s) of interest)";

  const char* gffFileFlag_ = "gffFile";
  const char* gffFileShortFlag_ = "g";

  addOption(parser, seqan::ArgParseOption(gffFileShortFlag_, gffFileFlag_, gff_desc, seqan::ArgParseArgument::INPUT_FILE, "GFF_FILE"));

  const char* fasta_desc =
  R"(The fasta file that contains the reference nucleotide sequence for the chromsome(s) of interest)";

  const char* fastaFileFlag_ = "fastaFile";
  const char* fastaFileShortFlag_ = "f";

  addOption(parser, seqan::ArgParseOption(fastaFileShortFlag_, fastaFileFlag_, fasta_desc, seqan::ArgParseArgument::INPUT_FILE, "FASTA_FILE"));

  const char* gaf_desc =
  R"(The gaf file that contains that annotates sequences defined in gff file. This is an optional parameter.)";

  const char* gafFileFlag_ = "gafFile";
  const char* gafFileShortFlag_ = "i";

  addOption(parser, seqan::ArgParseOption(gafFileShortFlag_, gafFileFlag_, gaf_desc, seqan::ArgParseArgument::INPUT_FILE, "GAF_FILE (OPT)"));


  const char* file_list_desc =
  R"(A file with a list of all SAM files to be read and analyzed. The SAM file reads should be
   aligned with respect to the reference genome specified by the --fasta flag.)";

  const char* fileListFlag_ = "fileList";
  const char* fileListShortFlag_ = "r";

  addOption(parser, seqan::ArgParseOption(fileListShortFlag_, fileListFlag_, file_list_desc, seqan::ArgParseArgument::INPUT_FILE, "SAM_FILE_LIST"));

  const char* contig_desc =
  R"(Define which contiguous DNA region (chromosome/mitochondria) to process.
  Defaults to '*' for all contiguous regions.)";

  const char* contigFlag_ = "contig";
  const char* contigShortFlag_ = "w";

  addOption(parser, seqan::ArgParseOption(contigShortFlag_, contigFlag_, contig_desc, seqan::ArgParseArgument::STRING, "CONTIG_ID"));


  const char* analysis_desc =
  R"(Define which genome analytic to perform. Defaults to '*' for all analytics.)";

  const char* analysisFlag_ = "analysisType";
  const char* analysisShortFlag_ = "z";

  addOption(parser, seqan::ArgParseOption(analysisShortFlag_, analysisFlag_, analysis_desc, seqan::ArgParseArgument::STRING, "ANALYSIS_TYPE"));


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

  const char* outCSV_desc =
  R"(Output all variant to an comma delimiter CSV file for further processing, (default "kgl_variant.csv").
  This file is always produced in the work directory.)";

  const char* CSVFileFlag_ = "outCSVFile";
  const char* CSVFileShortFlag_ = "o";

  addOption(parser, seqan::ArgParseOption(CSVFileShortFlag_, CSVFileFlag_, outCSV_desc, seqan::ArgParseArgument::OUTPUT_FILE, "OUT_CSV_FILE"));

  const char* min_count_desc =
  R"(The minimum SAM/BAM count of a single nucleotide analyzed in the Parent (wild-type)
  genome that is at variance from the reference (fasta) genome)";

  const char* minimumCountFlag_ = "minimumCount";
  const char* minimumCountShortFlag_ = "c";

  addOption(parser, seqan::ArgParseOption(minimumCountShortFlag_, minimumCountFlag_, min_count_desc , seqan::ArgParseArgument::INTEGER, "INT"));

  const char* min_prop_desc =
  R"(The minimum proportion of a single nucleotide analyzed in the Parent (wild-type)
  genome that is at variance from the reference (fasta) genome)";

  const char* minimumProportionFlag_ = "minimumProportion";
  const char* minimumProportionShortFlag_ = "p";

  addOption(parser, seqan::ArgParseOption(minimumProportionShortFlag_, minimumProportionFlag_, min_prop_desc, seqan::ArgParseArgument::DOUBLE, "FLOAT"));

  const char* translation_table_desc =
  R"(The amino acid translation table name used in gene verification and amino acid generation from CDS regions.
  The tables match the translation tables defines on the NCBI website at
  https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi. Individual organism tables are also supported.
  The current defined tables are: 'NCBI_TABLE_1' (the standard default table), 'NCBI_TABLE_2', 'NCBI_TABLE_3',
  'NCBI_TABLE_4', 'NCBI_TABLE_5' and organism 'P_FALCIPARUM')";

  const char* translationTableFlag_ = "translationTable";
  const char* translationTableShortFlag_ = "t";

  addOption(parser, seqan::ArgParseOption(translationTableShortFlag_, translationTableFlag_, translation_table_desc, seqan::ArgParseArgument::STRING, "STRING"));

  const char* gatk_desc =
  R"(Flag. ALL vcf files are parsed as formatted by the GATK format caller. Any FreeBayes generated vcf files will FAIL
  to parse if this flag is used. Note that if this flag is NOT specified, vcf file names tagged with the prefix
  'gatk<Name>.vcf' (case insensitive) are automatically treated as GATK formatted files.
  This allows a mixture of FreeBayes and GATK derived vcf files to be used in a single analysis.)";

  const char* gatkFlag_ = "gatk_only";
  const char* gatkShortFlag_ = "k";

  addOption(parser, seqan::ArgParseOption(gatkShortFlag_, gatkFlag_, gatk_desc, seqan::ArgParseArgument::BOOL, "FLAG"));


  const char* verbose_desc =
  R"(Flag. Enables verbose logged output to screen and log file.)";

  const char* verboseFlag_ = "verbose";
  const char* verboseShortFlag_ = "v";

  addOption(parser, seqan::ArgParseOption(verboseShortFlag_, verboseFlag_, verbose_desc, seqan::ArgParseArgument::BOOL, "FLAG"));

  const char* read_quality_desc =
  R"(The nucleotide read quality as -10 log10 Pr {ReadError} e.g. 30 is a 1/1000 chance
  of an nucleotide read error. Defaults to 30.)";

  const char* readQualityFlag_ = "readQuality";
  const char* readQualityShortFlag_ = "q";

  addOption(parser, seqan::ArgParseOption(readQualityShortFlag_, readQualityFlag_, read_quality_desc, seqan::ArgParseArgument::DOUBLE, "FLOAT"));

  const char* variant_quality_desc =
  R"(The variant quality as -10 log10 Pr {VariantError} e.g. 10 is a 1/10 chance the variant is incorrectly called. Defaults to 10.)";

  const char* variantQualityFlag_ = "variantQuality";
  const char* variantQualityShortFlag_ = "a";

  addOption(parser, seqan::ArgParseOption(variantQualityShortFlag_, variantQualityFlag_, variant_quality_desc, seqan::ArgParseArgument::DOUBLE, "FLOAT"));

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
  createLogger(MODULE_NAME, args().logFile, args().max_error_count, args().max_warn_count);

  std::string file_list;

  // Setup the files.
  getFilePath(fastaFileFlag_ ,parser , directory_path, args_.fastaFile);
  getFilePath(gffFileFlag_ ,parser , directory_path, args_.gffFile);

  if (seqan::isSet(parser, fileListFlag_)) {

    getFilePath(fileListFlag_ ,parser , directory_path, file_list);

  }

  if (file_list.length() == 0) {

    ExecEnv::log().critical("No SAM/BAM file list specified with --fileList");

  } else {

    if (not readFileList(file_list, directory_path, args_.fileList)) {

      ExecEnv::log().critical("Unable to read specified file list: {}", file_list);

    } else if (args_.fileList.empty()){

      ExecEnv::log().critical("No SAM/BAM files read from file list: {}", file_list);

    }

  }

  // Setup the gaf File (if specified).
  if (seqan::isSet(parser, gafFileFlag_)) {

    std::string gffFile;
    seqan::getOptionValue(gffFile, parser, gafFileFlag_);
    fs::path gff_file_path = directory_path / fs::path(gffFile);
    args_.gafFile = gff_file_path.string();

  }

  // Setup the outCSVFile.
  if (seqan::isSet(parser, CSVFileFlag_)) {

    std::string outCSVFile;
    seqan::getOptionValue(outCSVFile, parser, CSVFileFlag_);
    fs::path out_file_path = directory_path / fs::path(outCSVFile);
    args_.outCSVFile = out_file_path.string();

  } else {

    fs::path out_file_path = directory_path / fs::path(args_.outCSVFile);
    args_.outCSVFile = out_file_path.string();

  }
  // Join the log file and the directory
  // check by opening.
  std::fstream out_file(args_.outCSVFile, std::fstream::out | std::fstream::app);
  if (!out_file) {

    ExecEnv::log().critical("Cannot open output CSV file (--outCSVFile): {}", args_.outCSVFile);

  }

  // Setup Other Parameters.
  if (seqan::isSet(parser, contigFlag_)) seqan::getOptionValue(args_.contig, parser, contigFlag_);

  if (seqan::isSet(parser, analysisFlag_)) {

    seqan::getOptionValue(args_.analysisType, parser, analysisFlag_);

    args_.analysisType = Utility::toupper(args_.analysisType);

    if (args_.analysisType != Phylogenetic::WILDCARD
        and args_.analysisType != Phylogenetic::ANALYZE_INTERVAL
        and  args_.analysisType != Phylogenetic::ANALYZE_SEQUENCES
        and args_.analysisType != Phylogenetic::ANALYZE_GENE
        and args_.analysisType != Phylogenetic::ANALYZE_REGION
        and args_.analysisType != Phylogenetic::ANALYZE_UPGMA) {

      ExecEnv::log().critical("Invalid Analysis Type: {}.  Must be one of: {}, {}, {}, {}, {}.",
                              args_.analysisType,
                              Phylogenetic::ANALYZE_INTERVAL,
                              Phylogenetic::ANALYZE_SEQUENCES,
                              Phylogenetic::ANALYZE_GENE,
                              Phylogenetic::ANALYZE_REGION,
                              Phylogenetic::ANALYZE_UPGMA);

    }

  }

  if (seqan::isSet(parser, minimumCountFlag_)) seqan::getOptionValue(args_.minCount, parser, minimumCountFlag_);

  if (seqan::isSet(parser, minimumProportionFlag_)) seqan::getOptionValue(args_.minProportion, parser, minimumProportionFlag_);

  if (seqan::isSet(parser, translationTableFlag_)) seqan::getOptionValue(args_.aminoTranslationTable, parser, translationTableFlag_);

  if (seqan::isSet(parser, gatkFlag_)) seqan::getOptionValue(args_.vcfAllGATK, parser, gatkFlag_);

  if (seqan::isSet(parser, verboseFlag_)) seqan::getOptionValue(args_.verbose, parser, verboseFlag_);

  if (seqan::isSet(parser, readQualityFlag_)) seqan::getOptionValue(args_.readQuality, parser, readQualityFlag_);

  if (seqan::isSet(parser, variantQualityFlag_)) seqan::getOptionValue(args_.variantQuality, parser, variantQualityFlag_);

  ExecEnv::log().SetVerbose(args().verbose);  // Set the logger verbose level.

  return true;

}

