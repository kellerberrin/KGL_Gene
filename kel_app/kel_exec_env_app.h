//
// Created by kellerberrin on 22/10/20.
//

#ifndef KEL_EXEC_ENV_APP_H
#define KEL_EXEC_ENV_APP_H

#include "kel_exec_env.h"
#include "kel_utility.h"

#include <iostream>
#include <csignal>


namespace kellerberrin {   //  organization level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Specifies the application runtime object.
// This header file must be included in the source file containing the "main" program entry point.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// The runApplication() function sets up the static environment object.
/// and provides logging and commandline arguments to the
/// application object. The function is a template so that
/// different applications may be easily specified.

template<class Environment>
int ExecEnv::runApplication(int argc, char const ** argv) {


  try {

    // Save the command line.
    getCommandLine(argc, argv);

    // Setup the static ExecEnv runtime environment and create the logger.
    if (not Environment::parseCommandLine(argc, argv)) {

      std::cerr << Environment::MODULE_NAME << " " << Environment::VERSION << "ExecEnv::runApplication - cannot parse command line" << std::endl;
      std::exit(EXIT_FAILURE);

    }

    // Create the application logger
    log_ptr_ = Environment::createLogger();
    if (not log_ptr_) {

      std::cerr << Environment::MODULE_NAME << " " << Environment::VERSION << "ExecEnv::runApplication - cannot create application logger" << std::endl;
      std::exit(EXIT_FAILURE);

    }

    // Enable ctrl-C runtime termination.
    signal(SIGINT, ctrlC);

    log().info("############ {} {} Start Runtime ###########", Environment::MODULE_NAME, Environment::VERSION);
    log().info("Command Line: {}", commandLine());

    // Run the application.
    Environment::executeApp();

    //  Runtime stats.
    double Clock, System, User;
    Utility::getRuntime(Clock, System, User);
    log().info("Runtime seconds; Clock: {:.2f}, System CPU: {:.2f}, User CPU: {:.2f} (No GPU)", Clock, System, User);
    log().info("############ {} {} End Runtime ###########", Environment::MODULE_NAME, Environment::VERSION);

    // Explicitly shutdown the logger.
    log_ptr_ = nullptr;

  } catch(std::exception& e) { // In general, unhandled exceptions should not appear here, so complain and exit.

    std::cerr << Environment::MODULE_NAME << " " << Environment::VERSION << "ExecEnv::runApplication - Unexpected/Uncaught Exception: " << e.what() << std::endl;
    std::exit(EXIT_FAILURE);

  }

  return EXIT_SUCCESS;

}


} // namespace

#endif //KEL_EXEC_ENV_APP_H
