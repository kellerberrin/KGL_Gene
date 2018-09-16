//
// Created by kellerberrin on 4/05/18.
//

#ifndef KGL_PHYLOGENETIC_APP_H
#define KGL_PHYLOGENETIC_APP_H


#include "kgl_genome_types.h"
#include "kgl_exec_env.h"
#include "kgl_phylogenetic_option.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

struct FileListInfo {

  std::string file_name;
  std::string genome_name;

};

// Holds the Phylogenetic Arguments.
struct Phylogenetic {

  std::string workDirectory{"./Work"};
  std::string fastaFile{""};
  std::string gffFile{""};
  std::string gafFile{""};
  std::vector<FileListInfo> fileList;
  std::string auxCSVFile{"kgl_aux.csv"};
  std::string logFile{"kgl_phylo.log"};
  std::string contig{WILDCARD};
  std::string aminoTranslationTable{"NCBI_TABLE_1"};
  std::string analysisType{WILDCARD};
  size_t ploidy{1};
  int max_error_count{1000};
  int max_warn_count{1000};
  bool verbose{false};

  static constexpr const char* WILDCARD = "*";  // All analytics

};

// Singleton. The Phylogenetic Runtime environment.
class PhylogeneticExecEnv {

public:

  PhylogeneticExecEnv()=delete;
  ~PhylogeneticExecEnv()=delete;

// The following 4 static members are required for all applications.
  static constexpr const char* VERSION = "0.1";
  static constexpr const char* MODULE_NAME = "kgl_phylo";
  static void executeApp(); // Application mainline.
  static bool parseCommandLine(int argc, char const ** argv);  // Parse command line arguments.

private:

  static const Phylogenetic& getArgs();
  static const RuntimeOptions& getRuntimeOptions();
  static Phylogenetic args_;
  static RuntimeOptions runtime_options_;

};



} //  organization level namespace
}  // project level namespace



#endif //KGL_PHYLOGENETIC_APP_H
