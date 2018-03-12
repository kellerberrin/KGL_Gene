//
// Created by kellerberrin on 10/11/17.
//

#ifndef KGL_PHYLOGENETIC_ENV_H
#define KGL_PHYLOGENETIC_ENV_H


#include "kgl_genome_types.h"
#include "kgl_exec_env.h"

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
  int max_error_count{1000};
  int max_warn_count{1000};
  bool verbose{false};
  int minCount{0};
  double minProportion{0};
  Phred_t readQuality{30.0};
  Phred_t variantQuality{10.0};

  static constexpr const char* WILDCARD = "*";  // All analytics and All contigs.

  // Analytic types.
  static constexpr const char* ANALYZE_INTERVAL = "INTERVAL";
  static constexpr const char* ANALYZE_SEQUENCES = "SEQUENCE";
  static constexpr const char* ANALYZE_REGION = "REGION";
  static constexpr const char* ANALYZE_UPGMA = "UPGMA";
  static constexpr const char* ANALYZE_GENE = "GENE";
  static constexpr const char* ANALYZE_RNA = "RNA";

};

// Singleton. The Phylogenetic Runtime environment.
class PhylogeneticExecEnv : public ExecEnv {

public:

  PhylogeneticExecEnv()=delete;
  ~PhylogeneticExecEnv()=delete;

  static const Phylogenetic& args();
  static bool parseCommandLine(int argc, char const ** argv);
  static constexpr const char* VERSION = "0.1";
  static constexpr const char* MODULE_NAME = "kgl_phylo";


private:

  static Phylogenetic args_;

};



} //  organization level namespace
}  // project level namespace


#endif //KGL_PHYLOGENETIC_ENV_H
