//
// Created by kellerberrin on 30/09/17.
//
#include <csignal>
#include <iostream>
#include "kgl_utility.h"
#include "kgl_phylogenetic_env.h"
#include "kgl_phylogenetic_app.h"


namespace kgl = kellerberrin::genome;



void ctrlC(int) {

  kgl::ExecEnv::log().warn("Control-C. Program terminates. Output files may be corrupt. Multi-threaded code may hang.");
  std::exit(EXIT_FAILURE);

}

/// The application() function sets up the static environment object.
/// and provides logging and commandline arguments to the
/// application object. The function is a template so that
/// different applications may be easily specified.

template<class Environment, class Application>
int application(int argc, char const ** argv) {


  try {

    Environment::parseCommandLine(argc, argv);  // Setup the static ExecEnv runtime environment.

    signal(SIGINT, ctrlC);

    kgl::ExecEnv::log().info("############ {} {} Start Processing ###########", Environment::MODULE_NAME, Environment::VERSION);
    kgl::ExecEnv::log().info("Command Line: {}", kgl::ExecEnv::commandLine());


    Application(Environment::args()); // Run the application.

    double Clock, System, User;
    kgl::Utility::getElpasedTime(Clock, System, User);
    kgl::ExecEnv::log().info("Elapsed seconds; Clock: {}, System CPU: {}, User CPU: {} (No GPU)", Clock, System, User);
    kgl::ExecEnv::log().info("############ {} {} End Processing ###########", Environment::MODULE_NAME, Environment::VERSION);

  } catch(std::exception& e) { // Code should not throw any exceptions, so complain and exit.

    std::cerr << Environment::MODULE_NAME << " " << Environment::VERSION << " - Unexpected Exception: " << e.what() << std::endl;
    std::exit(EXIT_FAILURE);

  }

  return EXIT_SUCCESS;

}

/// The mainline.
int main(int argc, char const ** argv)
{

  return application<kgl::PhylogeneticExecEnv, kgl::PhylogeneticApp>(argc, argv);

}

