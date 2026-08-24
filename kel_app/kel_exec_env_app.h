//
// Copyright 2023 Kellerberrin
//

#ifndef KEL_EXEC_ENV_APP_H
#define KEL_EXEC_ENV_APP_H

#include "kel_exec_env.h"
#include "kel_utility.h"

#include <iostream>
#include <csignal>
#include <chrono>


namespace kellerberrin {   //  organization level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Specifies the application runtime object.
// This header file must be included in the source file containing the "main" program entry point.
//
// The runApplication() function sets up the static environment object.
// It provides logging and commandline arguments to the application object. The function is a template so that
// different applications may be easily specified.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/// runApplication() template. Sets up the static ExecEnv runtime environment, creates the logger,
/// parses the command line, and invokes Environment::executeApp().
///
/// The template parameter Environment must expose:
///   - a static parseCommandLine(int argc, char const**) -> bool
///   - a static createLogger() -> std::unique_ptr<ExecEnvLogger>
///   - a static executeApp()
///   - static constants MODULE_NAME and VERSION
template<class Environment>
int ExecEnv::runApplication(int argc, char const ** argv) {


  try {

    // Get the start time.
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
    // Save the command line.
    setCommandTokens(argc, argv);

    // Setup the static ExecEnv runtime environment and create the logger.
    if (not Environment::parseCommandLine(argc, argv)) {

      std::cerr << Environment::MODULE_NAME << " " << Environment::VERSION << " ExecEnv::runApplication - command line parse unsuccessful" << std::endl;
      std::exit(EXIT_FAILURE);

    }

    auto logger = Environment::createLogger();
    if (not logger) {

      std::cerr << Environment::MODULE_NAME << " " << Environment::VERSION << " ExecEnv::runApplication - cannot create application logger" << std::endl;
      std::exit(EXIT_FAILURE);

    }
    {
      std::lock_guard lock(env_mutex_);
      log_ptr_ = std::move(logger);
    }

    // Enable ctrl-C runtime termination.
    std::signal(SIGINT, ctrlC);

    log().info("############ {} {} Start Runtime ###########", Environment::MODULE_NAME, Environment::VERSION);
    log().info("Command Line: {}", commandLine());

    // Run the application.
    Environment::executeApp();

    //  Runtime stats.
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    double clock = static_cast<double>(elapsed.count()) / 1000.0;
    auto [system, user] = Utility::process_time_usage();
    log().info("Runtime seconds; Clock: {:.2f}, System CPU: {:.2f}, User CPU: {:.2f}", clock, system, user);
    log().info("############ {} {} End Runtime ###########", Environment::MODULE_NAME, Environment::VERSION);

    // Explicitly shutdown the logger.
    {
      std::lock_guard lock(env_mutex_);
      log_ptr_ = nullptr;
    }


  } catch(std::exception& e) { // In general, unhandled exceptions should not appear here, so complain and exit.

    std::cerr << Environment::MODULE_NAME << " " << Environment::VERSION << " ExecEnv::runApplication - Unexpected/Uncaught Exception: " << e.what() << std::endl;
    std::exit(EXIT_FAILURE);

  }

  return EXIT_SUCCESS;

}


} // namespace

#endif //KEL_EXEC_ENV_APP_H
