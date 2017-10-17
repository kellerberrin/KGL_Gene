// MIT License
//
// Copyright (c) 2017
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//
//
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
  std::string parentFile{"parent.sam"};
  std::string mutantFile{"mutant.sam"};
  std::string logFile{"kgl_snp.log"};
  std::string contig{"*"};
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
  static constexpr const char* MODULE_NAME = "kgl_minority";

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
