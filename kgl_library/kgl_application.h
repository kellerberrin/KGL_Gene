//
// Created by kellerberrin on 4/05/18.
//

#ifndef KGL_APPLICATION_H
#define KGL_APPLICATION_H

#include <csignal>
#include <iostream>
#include "kgl_utility.h"
#include "kgl_exec_env.h"



namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



void ctrlC(int) {

  ExecEnv::log().warn("Control-C. Program terminates. Output files may be corrupt. Multi-threaded code may hang.");
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

    ExecEnv::log().info("############ {} {} Start Processing ###########", Environment::MODULE_NAME, Environment::VERSION);
    ExecEnv::log().info("Command Line: {}", ExecEnv::commandLine());


    Application(Environment::args()); // Run the application.

    double Clock, System, User;
    Utility::getElpasedTime(Clock, System, User);
    ExecEnv::log().info("Elapsed seconds; Clock: {}, System CPU: {}, User CPU: {} (No GPU)", Clock, System, User);
    ExecEnv::log().info("############ {} {} End Processing ###########", Environment::MODULE_NAME, Environment::VERSION);

  } catch(std::exception& e) { // Code should not throw any exceptions, so complain and exit.

    std::cerr << Environment::MODULE_NAME << " " << Environment::VERSION << " - Unexpected Exception: " << e.what() << std::endl;
    std::exit(EXIT_FAILURE);

  }

  return EXIT_SUCCESS;

}



}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_APPLICATION_H
