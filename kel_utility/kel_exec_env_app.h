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

    // Create the template execution environment
    auto environment_ptr = std::make_unique<Environment>();

    // Setup the static ExecEnv runtime environment and create the logger.
    if (not environment_ptr->parseCommandLine(argc, argv)) {

      std::cerr << Environment::MODULE_NAME << " " << Environment::VERSION << "ExecEnv::runApplication - cannot parse command line" << std::endl;
      std::exit(EXIT_FAILURE);

    }

    signal(SIGINT, ctrlC);

    log().info("############ {} {} Start Processing ###########", Environment::MODULE_NAME, Environment::VERSION);
    log().info("Command Line: {}", commandLine());

    environment_ptr->executeApp(); // Run the application.

    double Clock, System, User;
    Utility::getElapsedTime(Clock, System, User);
    log().info("Elapsed seconds; Clock: {:.2f}, System CPU: {:.2f}, User CPU: {:.2f} (No GPU)", Clock, System, User);
    log().info("############ {} {} End Processing ###########", Environment::MODULE_NAME, Environment::VERSION);

    environment_ptr = nullptr; // shutdown the application
    log_ptr_ = nullptr; // shutdown the logger.

  } catch(std::exception& e) { // Code should not throw any unhandled exceptions, so complain and exit.

    std::cerr << Environment::MODULE_NAME << " " << Environment::VERSION << "ExecEnv::runApplication - Unexpected/Uncaught Exception: " << e.what() << std::endl;
    std::exit(EXIT_FAILURE);

  }

  return EXIT_SUCCESS;

}


} // namespace

#endif //KEL_EXEC_ENV_APP_H
