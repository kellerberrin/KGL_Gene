// MIT License
//
// Copyright (c) 2017
//
// Created by kellerberrin on 17/10/17.
//

#ifndef KGL_MINORITY_ENV_H
#define KGL_MINORITY_ENV_H

#include "kgl_exec_env.h"

namespace kellerberrin {   //  organization level namespace
  namespace genome {   // project level namespace


// Holds the Minority Arguments.
struct MinorityArgs {

  std::string workDirectory{"./Work"};
  std::string fastaFile{""};
  std::string gffFile{""};
  std::string parentFile{""};
  std::vector<std::string> fileList;
  std::string mutantFile{""};
  std::string logFile{"kgl_snp.log"};
  std::string contig{"*"};
  size_t aminoTranslationTable{1};
  bool verbose{false};
  int mutantMinCount{20};
  double mutantMinProportion{0.7};
  int parentMinCount{20};
  double parentMinProportion{0.7};
  int threadCount{-1};
  unsigned char readQuality{0};
  int lockGranularity{1000};
  int queueSize{1000000};

};

// The Minority Runtime environment.
class MinorityExecEnv : public ExecEnv {

public:

  MinorityExecEnv()=delete;
  ~MinorityExecEnv()=delete;

  static const MinorityArgs& args();
  static bool parseCommandLine(int argc, char const ** argv);
  static constexpr const char* VERSION = "0.1";
  static constexpr const char* MODULE_NAME = "kgl_phylo";

  class Application;

private:

  static MinorityArgs args_;

};

// Application class implements the mainline logic and controls
// data object lifetimes, see kgl_minority_app.cc.
class MinorityExecEnv::Application {

public:

  Application(Logger& log, const MinorityArgs& args);
  ~Application() = default;

};



  }; //  organization level namespace
};  // project level namespace

#endif //KGL_MINORITY_ENV_H
