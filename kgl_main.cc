//
// Created by kellerberrin on 30/09/17.
//
#include <csignal>
#include <iostream>
#include "kgl_phylogenetic_env.h"

namespace kgl = kellerberrin::genome;


void ctrlC(int sig) {

  kgl::ExecEnv::log().warn("Control-C. Program terminates. Output files may be corrupt. Multi-threaded code may hang.");
  std::exit(EXIT_FAILURE);

}


// The application() function sets up the static environment object.
// and provides logging and commandline arguments to the
// application object. The function is a template so that
// different applications may be easily specified.

template<class T>
int application(int argc, char const ** argv) {


  try {

    T::parseCommandLine(argc, argv);  // Setup the static ExecEnv runtime environment.

    signal(SIGINT, ctrlC);

    kgl::ExecEnv::log().info("############ {} {} Start Processing ###########", T::MODULE_NAME, T::VERSION);

    typename T::Application(T::log(), T::args()); // Run the application.

    double Clock, System, User;
    kgl::ExecEnv::getElpasedTime(Clock, System, User);
    kgl::ExecEnv::log().info("Elapsed seconds; Clock: {}, System CPU: {}, User CPU: {} (No GPU)", Clock, System, User);
    kgl::ExecEnv::log().info("############ {} {} End Processing ###########", T::MODULE_NAME, T::VERSION);

  } catch(std::exception& e) { // Code should not throw any exceptions, so complain and exit.

    std::cerr << T::MODULE_NAME << " " << T::VERSION << " - Unexpected Exception: " << e.what() << std::endl;
    std::exit(EXIT_FAILURE);

  }

  return EXIT_SUCCESS;

}


int main(int argc, char const ** argv)
{

  return application<kgl::PhylogeneticExecEnv>(argc, argv);

}

