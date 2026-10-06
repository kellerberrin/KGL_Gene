//
// Created by kellerberrin on 4/05/18.
//

#ifndef KGL_GENE_APP_H
#define KGL_GENE_APP_H


#include "kel_exec_env.h"


namespace kellerberrin::genome {   //  organization::project level namespace


/// Holds the Commandline Arguments, initialized with defaults.
struct CmdLineArgs {

  /// Working directory for file paths.
  std::string reference_directory{"./"};
  /// Log file name.
  std::string log_file{"kgl_phylo.log"};
  /// Runtime options XML file name.
  std::string options_file{"runtime_options.xml"};
  /// Optional file for the post-processed options_file.
  std::string option_file_out;
  /// Optional file output of parsed XML options.
  std::string parsed_option_out;
  /// Maximum number of error messages before suppression.
  size_t max_error_count{1000};
  /// Maximum number of warning messages before suppression.
  size_t max_warn_count{1000};

};

/// The Runtime environment.
class GeneExecEnv {

public:

  GeneExecEnv()=delete;
  ~GeneExecEnv()=delete;


// The following 5 static members are required for all applications.
  /// Application version string.
  inline static constexpr const char* VERSION = "1.0";
  /// Module name string.
  inline static constexpr const char* MODULE_NAME = "kglGene";
  /// Application mainline.
  static void executeApp();
  /// Parse command line arguments.
  [[nodiscard]] static bool parseCommandLine(int argc, char const ** argv);
  /// Create application logger.
  [[nodiscard]] static std::unique_ptr<ExecEnvLogger> createLogger();
  /// Returns the command line arguments.
  [[nodiscard]] static const CmdLineArgs& getArgs() { return args_; }

private:

  inline static CmdLineArgs args_;

};



} //  end namespace



#endif //KGL_GENE_APP_H