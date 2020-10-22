///
// Created by kellerberrin on 30/09/17.
//

#include "kel_exec_env.h"
#include "kel_utility.h"


// Define namespace alias
namespace kel = kellerberrin;



class ImplLogger {

public:

  ImplLogger( const std::string& module,
              const std::string& log_file,
              int max_error_messages,
              int max_warning_messages);
  ~ImplLogger() = default;

  template<typename... Args> void trace(const char* message, Args... args) noexcept { log_ptr_->trace(message, args...); }
  template<typename... Args> void info(const char* message, Args... args) noexcept { log_ptr_->info(message, args...); }
  template<typename... Args> void warn(const char* message, Args... args) noexcept { log_ptr_->warn(message, args...); }
  template<typename... Args> void error(const char* message, Args... args) noexcept { log_ptr_->error(message, args...); }
  template<typename... Args> void critical(const char* message, Args... args) noexcept { log_ptr_->critical(message, args...); }
  // These function only report if the verbose flag is set; see "setVerbose(true)" above.
  template<typename... Args> void vtrace(const char* message, Args... args) noexcept { log_ptr_->vtrace(message, args...); }
  template<typename... Args> void vinfo(const char* message, Args... args) noexcept { log_ptr_->vinfo(message, args...); }
  template<typename... Args> void vwarn(const char* message, Args... args) noexcept { log_ptr_->vwarn(message, args...); }

private:

  std::unique_ptr<kel::Logger> log_ptr_;

};



ImplLogger::ImplLogger( const std::string& module,
                        const std::string& log_file,
                        int max_error_messages,
                        int max_warning_messages) {

  log_ptr_ = std::make_unique<kel::Logger>(module, log_file);

  log_ptr_->SetMaxErrorMessages(max_error_messages);

  log_ptr_->SetMaxwarningMessages(max_warning_messages);

}


// The ExecEnv private static variables.
std::string kel::ExecEnv::command_line_;
std::unique_ptr<kel::Logger> kel::ExecEnv::log_ptr_;
std::unique_ptr<ImplLogger> kel::ExecEnv::impl_log_ptr_;


void kel::ExecEnv::createLogger(const std::string& module,
                                const std::string& log_file,
                                int max_error_messages,
                                int max_warning_messages) {

  log_ptr_ = std::make_unique<Logger>(module, log_file);

  log_ptr_->SetMaxErrorMessages(max_error_messages);

  log_ptr_->SetMaxwarningMessages(max_warning_messages);

//  impl_log_ptr_ = std::make_unique<ImplLogger>(module, log_file, max_error_messages, max_warning_messages);

}


void kel::ExecEnv::ctrlC(int) {

  ExecEnv::log().warn("Control-C. Program terminates. Output files may be corrupt. Multi-threaded code may hang.");
  std::exit(EXIT_FAILURE);

}

void kel::ExecEnv::getCommandLine(int argc, char const ** argv) {

  std::string command_line;

  for (int idx = 0; idx < argc; ++idx) {

    command_line += argv[idx];
    command_line += ' ';

  }

  command_line_ = command_line;

}


