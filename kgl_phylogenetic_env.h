//
// Created by kellerberrin on 10/11/17.
//

#ifndef KGL_PHYLOGENETIC_ENV_H
#define KGL_PHYLOGENETIC_ENV_H


#include "kgl_exec_env.h"

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


// Holds the Minority Arguments.
struct Phylogenetic {

  std::string workDirectory{"./Work"};
  std::string fastaFile{""};
  std::string gffFile{""};
  std::vector<std::string> fileList;
  std::string outCSVFile{"kgl_out.csv"};
  std::string logFile{"kgl_snp.log"};
  std::string contig{"*"};
  size_t aminoTranslationTable{1};
  int max_error_count{1000};
  int max_warn_count{1000};
  bool verbose{false};
  int minCount{0};
  double minProportion{0};
  unsigned char readQuality{0};

};

// The Phylogenetic Runtime environment.
class PhylogeneticExecEnv : public ExecEnv {

public:

  PhylogeneticExecEnv()=delete;
  ~PhylogeneticExecEnv()=delete;

  static const Phylogenetic& args();
  static bool parseCommandLine(int argc, char const ** argv);
  static constexpr const char* VERSION = "0.1";
  static constexpr const char* MODULE_NAME = "kgl_phylo";

  class Application;

private:

  static Phylogenetic args_;

};

// Application class implements the mainline logic and controls
// data object lifetimes, see kgl_phylogenetic_app.cc.
class PhylogeneticExecEnv::Application {

public:

  Application(Logger& log, const Phylogenetic& args);
  ~Application() = default;

};



}; //  organization level namespace
};  // project level namespace


#endif //KGL_PHYLOGENETIC_ENV_H
