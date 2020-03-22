//
// Created by kellerberrin on 30/09/17.
//

#ifndef KGL_EXEC_ENV_H
#define KGL_EXEC_ENV_H

#include <csignal>
#include <iostream>
#include <string>
#include <memory>
#include "kel_logging.h"
#include "kel_utility.h"


namespace kellerberrin {   //  organization level namespace

// Singleton. This class sets up the application runtime environment as a series of static variables
// and member functions. The class is never instantiated and is the first and only statement in main() (see kgl_main.cc).

class ExecEnv {

public:

  ExecEnv()=delete;
  ~ExecEnv()=delete;


  static Logger& log();
  template<class Environment> static int runApplication(int argc, char const ** argv);

  static void getCommandLine(int argc, char const ** argv);
  static void createLogger(const std::string& module,
                           const std::string& log_file,
                           int max_error_message,
                           int max_warning_messages);
  static const std::string& commandLine() { return command_line_; }


private:

  static std::string command_line_;
  static std::unique_ptr<Logger> log_ptr_;

  static void ctrlC(int);

};



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
    Environment environment;

    environment.parseCommandLine(argc, argv);  // Setup the static ExecEnv runtime environment and create the logger.

    signal(SIGINT, ctrlC);

    log().info("############ {} {} Start Processing ###########", Environment::MODULE_NAME, Environment::VERSION);
    log().info("Command Line: {}", commandLine());


    environment.executeApp(); // Run the application.

    double Clock, System, User;
    Utility::getElapsedTime(Clock, System, User);
    log().info("Elapsed seconds; Clock: {}, System CPU: {}, User CPU: {} (No GPU)", Clock, System, User);
    log().info("############ {} {} End Processing ###########", Environment::MODULE_NAME, Environment::VERSION);

  } catch(std::exception& e) { // Code should not throw any exceptions, so complain and exit.

    std::cerr << Environment::MODULE_NAME << " " << Environment::VERSION << "ExecEnv::runApplication - Unexpected/Uncaught Exception: " << e.what() << std::endl;
    std::exit(EXIT_FAILURE);

  }

  return EXIT_SUCCESS;

}


}   // end namespace

#endif //KGL_EXEC_ENV_H
