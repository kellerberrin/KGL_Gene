//
// Created by kellerberrin on 8/05/18.
//


#include "kgd_deploid_app.h"
#include "kgl_utility.h"

#include <iostream>
#include <ctime>
#include <cctype>
#include <seqan/arg_parse.h>
#define BOOST_FILESYSTEM_NO_DEPRECATED // Recommended by boost filesystem documentation.
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>


namespace kgl = kellerberrin::genome;
namespace kgd = kellerberrin::deploid;
namespace fs = boost::filesystem;
namespace bt = boost;


// Utility function to pre-pend work directory path and check file existance
void getFilePath(const std::string& option_text,
                 seqan::ArgumentParser& parser,
                 fs::path& directory_path,
                 std::string& file_name,
                 bool isInputFile) {

  if (seqan::isSet(parser, option_text)) {

    seqan::getOptionValue(file_name, parser, option_text);

  } else if (file_name == kgd::DeploidArgs::NOT_SPECIFIED) {

    return;  // Optional argument so just return.

  }

  fs::path file_path = directory_path / fs::path(file_name);

  file_name = file_path.string();

  if (not isInputFile) return;

  // If an input file then check for existance.
  boost::system::error_code error_code;
  bool file_exists = fs::exists(file_path, error_code);

  if (!file_exists) {

    kgl::ExecEnv::log().critical("Option --{}, File: {} does not exist.", option_text, file_path.string());

  }

  if (error_code.value() != boost::system::errc::success) {

    kgl::ExecEnv::log().critical("Option --{}, File: {}, error: {}.", option_text, file_path.string(),  error_code.message());

  }


}


bool kgd::DeploidExecEnv::parseCommandLine(int argc, char const ** argv) {

  // Save the command line.
  kgl::ExecEnv::getCommandLine(argc, argv);

  // Setup ArgumentParser.
  seqan::ArgumentParser parser(MODULE_NAME);
  // Set short description, version, and date.
  setShortDescription(parser, "kgd_deploid COI haplotype generation");
  setVersion(parser, VERSION);
  setDate(parser, "May 2018");

  const char* program_desc =
  R"(kgd_deploid Complexity Of Infection (COI) haplotype generation)";

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

  addOption(parser, seqan::ArgParseOption(plafFileShortFlag_, plafFileFlag_, plaf_file_desc, seqan::ArgParseArgument::INPUT_FILE, "PLAF_FILE"));

  const char* outfile_file_desc =
  R"(The output file spec, without file extension, used as a template for multiple output files.
This is relative to the work directory the workDirectory is prepended to to this path e.g. 'workDirectory/output_file_path<extension>')";

  const char* outFileFlag_ = "outFile";
  const char* outFileShortFlag_ = "o";

  addOption(parser, seqan::ArgParseOption(outFileShortFlag_, outFileFlag_, outfile_file_desc, seqan::ArgParseArgument::STRING, "OUTPUT"));

  const char* panel_file_desc =
  R"(Strain panel file. This is relative to the work directory. and the workDirectory
  is prepended to to this path e.g. 'workDirectory/strain_panel_file_path')";

  const char* panelFileFlag_ = "strainPanelFile";
  const char* panelFileShortFlag_ = "panel";

  addOption(parser, seqan::ArgParseOption(panelFileShortFlag_, panelFileFlag_, panel_file_desc, seqan::ArgParseArgument::INPUT_FILE, "PANEL_FILE"));

  const char* exclude_allele_file_desc =
  R"(Excluded alleles file. This is relative to the work directory. and the workDirectory
  is prepended to to this path e.g. 'workDirectory/exclude_allele_file_path')";

  const char* excludeFileFlag_ = "excludedAllelesFile";
  const char* excludeFileShortFlag_ = "exclude";

  addOption(parser, seqan::ArgParseOption(excludeFileShortFlag_, excludeFileFlag_, exclude_allele_file_desc, seqan::ArgParseArgument::INPUT_FILE, "EXCLUDE_FILE"));


  const char* paint_strain_file_desc =
  R"(Paint the calculated strains to file. This is relative to the work directory. and the workDirectory
  is prepended to to this path e.g. 'workDirectory/paint_strain_file_path')";

  const char* paintFileFlag_ = "paintStrainsFile";
  const char* paintFileShortFlag_ = "paint";

  addOption(parser, seqan::ArgParseOption(paintFileShortFlag_, paintFileFlag_, paint_strain_file_desc, seqan::ArgParseArgument::OUTPUT_FILE, "PAINT_FILE"));

  const char* vcf_out_file_desc =
  R"(Save final halpotypes into an output VCF file. This is relative to the work directory. and the workDirectory
  is prepended to to this path e.g. 'workDirectory/paint_strain_file_path')";

  const char* vcfOutFileFlag_ = "vcfOutputFile";
  const char* vcfOutShortFlag_ = "vcfOut";

  addOption(parser, seqan::ArgParseOption(vcfOutShortFlag_, vcfOutFileFlag_, vcf_out_file_desc, seqan::ArgParseArgument::OUTPUT_FILE, "VCF_OUT_FILE"));

  const char* strain_count_desc =
  R"(The maximun number of canonical strains considered. The algorithm runtime is proportion to the square of this parameter.)";

  const char* strainCountFlag_ = "maximumStrainCount";
  const char* strainCountShortFlag_ = "k";

  addOption(parser, seqan::ArgParseOption(strainCountShortFlag_, strainCountFlag_, strain_count_desc , seqan::ArgParseArgument::INTEGER, "INT"));

  const char* MCMC_sample_desc =
  R"(The number of MCMC samples (default is 800))";

  const char* MCMCSampleFlag_ = "MCMCSampleCount";
  const char* MCMCSampleShortFlag_ = "nSample";

  addOption(parser, seqan::ArgParseOption(MCMCSampleShortFlag_, MCMCSampleFlag_, MCMC_sample_desc , seqan::ArgParseArgument::INTEGER, "INT"));

  const char* MCMC_rate_desc =
  R"(The MCMC sample rate (default value 5))";

  const char* MCMCRateFlag_ = "MCMCSampleRate";
  const char* MCMCRateShortFlag_ = "rate";

  addOption(parser, seqan::ArgParseOption(MCMCRateShortFlag_, MCMCRateFlag_, MCMC_rate_desc , seqan::ArgParseArgument::INTEGER, "INT"));

  const char* MCMC_random_seed =
  R"(The MCMC kgd_random seed (default uses time()))";

  const char* MCMCRandomSeedFlag_ = "MCMCRandomSeed";
  const char* MCMCRandomSeedShortFlag_ = "seed";

  addOption(parser, seqan::ArgParseOption(MCMCRandomSeedShortFlag_, MCMCRandomSeedFlag_, MCMC_random_seed , seqan::ArgParseArgument::INTEGER, "INT"));

  const char* MCMC_burn_desc =
  R"(MCMC burn rate (default value 0.5))";

  const char* MCMCBurnFlag_ = "MCMCBurnRate";
  const char* MCMCBurnShortFlag_ = "burn";

  addOption(parser, seqan::ArgParseOption(MCMCBurnShortFlag_, MCMCBurnFlag_, MCMC_burn_desc , seqan::ArgParseArgument::DOUBLE, "FLOAT"));


  const char* nopanel_desc =
  R"(Flag. Generate the candidate strains using IDB functionality. Do not use a pre-determined panel of strains)";

  const char* nopanelFlag_ = "noPanel";
  const char* nopanelShortFlag_ = "noPanel";

  addOption(parser, seqan::ArgParseOption(nopanelShortFlag_, nopanelFlag_, nopanel_desc));

  const char* idb_desc =
  R"(Flag. Calculate the canonical strains using identity by descent.)";

  const char* ibdFlag_ = "identityByDescent";
  const char* ibdShortFlag_ = "ibd";

  addOption(parser, seqan::ArgParseOption(ibdShortFlag_, ibdFlag_, idb_desc));


  const char* idb_painting_desc =
  R"(Flag. IBD painting, compute posterior probabilities of IBD configurations of given strain proportions.
This option must be used with option -initialP.)";

  const char* ibdPaintingFlag_ = "identityByDescentPainting";
  const char* ibdPaintingShortFlag_ = "ibdPainting";

  addOption(parser, seqan::ArgParseOption(ibdPaintingShortFlag_, ibdPaintingFlag_, idb_painting_desc));

  const char* inbreeding_desc =
  R"(Flag. Calculate the inbreeding probabilities.)";

  const char* inBreedingFlag_ = "inbreedingProbabilities";
  const char* inBreedingShortFlag_ = "inbreeding";

  addOption(parser, seqan::ArgParseOption(inBreedingShortFlag_, inBreedingFlag_, inbreeding_desc));

  const char* strain_proportion_desc =
  R"(Initialize Strain proportions as a list of floats (e.g. --initialStrainProportions 0.2 0.8).)";

  const char* StrainProportionList_ = "initialStrainProportions";
  const char* StrainProportionListShort_ = "initialP";

  addOption(parser, seqan::ArgParseOption(StrainProportionList_, StrainProportionListShort_, strain_proportion_desc , seqan::ArgParseArgument::DOUBLE, "FLOAT", true, 2));

  ///////////////////////////////////////////////////////////////////
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
  kgl::ExecEnv::createLogger(MODULE_NAME, getArgs().logFile, getArgs().max_error_count, getArgs().max_warn_count);

  // Setup input files.
  getFilePath(vcfFileFlag_, parser, directory_path, args_.vcfFile, true);
  getFilePath(plafFileFlag_, parser, directory_path, args_.plafFile, true);
  getFilePath(panelFileFlag_, parser, directory_path, args_.panelFile, true);
  getFilePath(excludeFileFlag_, parser, directory_path, args_.excludeFile, true);
  // Setup output files.
  getFilePath(paintFileFlag_, parser, directory_path, args_.paintFile, false);
  getFilePath(outFileFlag_, parser, directory_path, args_.outputTemplate, false);
  getFilePath(vcfOutFileFlag_, parser, directory_path, args_.vcfOutFile, false);

  // Setup Other Parameters.
  if (seqan::isSet(parser, strainCountFlag_)) seqan::getOptionValue(args_.maxStrains, parser, strainCountFlag_);
  if (seqan::isSet(parser, MCMCSampleFlag_)) seqan::getOptionValue(args_.MCMCSamples, parser, MCMCSampleFlag_);
  if (seqan::isSet(parser, MCMCRateFlag_)) seqan::getOptionValue(args_.MCMCSampleRate, parser, MCMCRateFlag_);
  if (seqan::isSet(parser, MCMCBurnFlag_)) seqan::getOptionValue(args_.MCMCBurnRate, parser, MCMCBurnFlag_);
  if (seqan::isSet(parser, nopanelFlag_)) seqan::getOptionValue(args_.noPanelFlag, parser, nopanelFlag_);
  if (seqan::isSet(parser, ibdFlag_)) seqan::getOptionValue(args_.identityByDescentFlag, parser, ibdFlag_);
  if (seqan::isSet(parser, inBreedingFlag_)) seqan::getOptionValue(args_.inbreedingProbabilitiesFlag, parser, inBreedingFlag_);
  if (seqan::isSet(parser, ibdPaintingFlag_)) seqan::getOptionValue(args_.identityByDescentPainting, parser, ibdPaintingFlag_);
  if (seqan::isSet(parser, inBreedingFlag_)) seqan::getOptionValue(args_.inbreedingProbabilitiesFlag, parser, inBreedingFlag_);
  if (seqan::isSet(parser, StrainProportionList_)) seqan::getOptionValue(args_.initialStrainProportions, parser, StrainProportionList_);

  // noPanel and Panel are mutually exclusive
  if (seqan::isSet(parser, nopanelFlag_) and seqan::isSet(parser, panelFileFlag_)) {

    kgl::ExecEnv::log().critical("Command line error. Option --{} and Option --{} are mutually exclusive", nopanelFlag_,  panelFileFlag_);

  }

  // If seed is not specified then use time() to seed rng.
  if (seqan::isSet(parser, MCMCRandomSeedFlag_)) {

    seqan::getOptionValue(args_.MCMCRandomSeed, parser, MCMCRandomSeedFlag_);

  } else {

    args_.MCMCRandomSeed = std::time(nullptr);

  }

  if (seqan::isSet(parser, StrainProportionList_) and seqan::isSet(parser, strainCountFlag_)) {

    kgl::ExecEnv::log().critical("Command line error. Option --{} and Option --{} are mutually exclusive", StrainProportionList_,  strainCountFlag_);

  }

  return true;

}