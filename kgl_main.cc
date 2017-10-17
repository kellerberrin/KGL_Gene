//
// Created by kellerberrin on 30/09/17.
//
#include <csignal>
#include <iostream>
#include "kgl_minority_env.h"

namespace kgl = kellerberrin::genome;


void ctrlC(int sig) {

  kgl::ExecEnv::log().warn("Control-C pressed. Program terminates. Open files may be in an unsafe state.");
  std::exit(EXIT_FAILURE);

}


// The application() function sets up the static environment object.
// and provides logging and commandline arguments to the
// application object. The function is a template so that
// different applications may be easily specified.

template<class T>
int application(int argc, char const ** argv) {


  try {

    T::parseCommandLine(argc, argv);  // Setup the runtime environment.

    signal(SIGINT, ctrlC);

    kgl::ExecEnv::log().info("############ {} {} Start Processing ###########",
                             T::MODULE_NAME,
                             T::VERSION);

    typename T::Application(T::log(), T::args()); // Do the analysis.

    double Clock, System, User;
    kgl::ExecEnv::getElpasedTime(Clock, System, User);
    kgl::ExecEnv::log().info("Elapsed seconds; Clock: {}, System CPU: {}, User CPU: {} (No GPU)", Clock, System, User);
    kgl::ExecEnv::log().info("############ {} {} End Processing ###########",
                             T::MODULE_NAME,
                             T::VERSION);

  } catch(...) { // Code should not throw any exceptions, so complain and exit.

    std::cerr << T::MODULE_NAME << " " << T::VERSION << " - Unexpected Exception." << std::endl;
    std::exit(EXIT_FAILURE);

  }

  return EXIT_SUCCESS;

}


int main(int argc, char const ** argv)
{

  return application<kgl::MinorityExecEnv>(argc, argv);

}

