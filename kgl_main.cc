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
#include <signal.h>
#include <iostream>
#include "kgl_exec_env.h"
#include "kgl_genome.h"

namespace kgl = kellerberrin::genome;


void ctrlC(int sig) {

  kgl::ExecEnv::log().warn("Control-C pressed. Program terminates. Open files may be in an unsafe state.");
  std::exit(EXIT_FAILURE);

}


int main(int argc, char const ** argv)
{

  try {

    kgl::ExecEnv::parseCommandLine(argc, argv);  // Setup the runtime environment.

    signal(SIGINT, ctrlC);

    kgl::ExecEnv::log().info("############ {} {} Start Processing ###########",
                             kgl::ExecEnv::MODULE_NAME,
                             kgl::ExecEnv::VERSION);

    kgl::GenomeAnalysis(kgl::ExecEnv::log(), kgl::ExecEnv::args()); // Do the analysis.

    double Clock, System, User;
    kgl::ExecEnv::getElpasedTime(Clock,System,User);
    kgl::ExecEnv::log().info("Elapsed seconds; Clock: {}, System CPU: {}, User CPU: {} (No GPU)", Clock, System, User);
    kgl::ExecEnv::log().info("############ {} {} End Processing ###########",
                             kgl::ExecEnv::MODULE_NAME,
                             kgl::ExecEnv::VERSION);

  } catch(...) { // Code should not throw any exceptions, so complain and exit.

    std::cerr << kgl::ExecEnv::MODULE_NAME << " " << kgl::ExecEnv::VERSION << " - Unexpected Exception." << std::endl;
    std::exit(EXIT_FAILURE);

  }

  return EXIT_SUCCESS;
}

