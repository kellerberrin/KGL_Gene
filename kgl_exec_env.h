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
// Created by kellerberrin on 30/09/17.
//

#ifndef KGL_EXEC_ENV_H
#define KGL_EXEC_ENV_H

#include "kgl_logging.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

// This class sets up the application runtime environment as a series of static variables and member functions.
// The class is never instantiated and is the first statement in main() (see kgl_main.cc).
class ExecEnv {

public:

  ExecEnv()=delete;
  ~ExecEnv()=delete;

  struct Args {

    std::string workDirectory{"./Work"};
    std::string fastaFile{""};
    std::string gffFile{""};
    std::string parentFile{""};
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

  static bool parseCommandLine(int argc, char const ** argv);
  static const Args& args() { return args_; }
  static Logger& log() { return *log_ptr_; }
  static constexpr const char* VERSION = "0.1";
  static constexpr const char* MODULE_NAME = "kgl_genome";
  static void getElpasedTime(double& Clock, double& System, double& User);

private:

  static Args args_;
  static std::unique_ptr<Logger> log_ptr_;

};


}   // namespace genome
}   // namespace kellerberrin

#endif //KGL_EXEC_ENV_H
