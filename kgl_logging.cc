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
// Created by kellerberrin on 6/09/17.
//


#include <iostream>

#include "spdlog/spdlog.h"  // Implement the logger using the spdlog library
#include "kgl_logging.h"


namespace kgl = kellerberrin::genome;
namespace spd = spdlog;

kgl::Logger& kgl::log(void) { // Global function returns a static instance of the logger.

  static kgl::Logger logger;
  return logger;

}

kgl::Logger::Logger(const std::string& log_file) : plog_impl_{std::make_unique<LogImpl>()} {}

class kgl::Logger::LogImpl {

public:

  LogImpl(const std::string& log_file) {

    std::vector<spdlog::sink_ptr> sinks;
    sinks.push_back( std::make_shared<spdlog::sinks::stdout_sink_mt>() );
    sinks.push_back( std::make_shared<spdlog::sinks::simple_file_sink_mt>(log_file) );
    plogger_ = std::make_unique<spd::logger>("readsam",  sinks.begin(), sinks.end());
    plogger_->set_level( spd::level::info);

  }
  ~LogImpl() = default;

  friend class kgl::Logger;

private:

  std::shared_ptr<spd::logger> plogger_;

};

void kgl::Logger::SetLevel(Severity level) noexcept {

  switch(level) {

    case Severity::Trace:
      spd::set_level(spd::level::trace);
      break;

    case Severity::Info:
      spd::set_level(spd::level::info);
      break;

    case Severity::Warn:
      spd::set_level(spd::level::warn);
      break;

    case Severity::Error:
      spd::set_level(spd::level::err);
      break;

    case Severity::Critical:
      spd::set_level(spd::level::critical);
      break;

  }

}

void kgl::Logger::SetFormat(const std::string& format) noexcept {

  spd::set_pattern(format);

}

template<typename M, typename... Args> void kgl::Logger::trace(M& message, Args... args) noexcept {

  plog_impl_->plogger_->trace(message, args);

}

template<typename M, typename... Args> void kgl::Logger::info(M& message, Args... args) noexcept {

  plog_impl_->plogger_->info(message, args);

}

template<typename M, typename... Args> void kgl::Logger::warn(M& message, Args... args) noexcept {

  plog_impl_->plogger_->warn(message, args);

}

template<typename M, typename... Args> void kgl::Logger::error(M& message, Args... args) noexcept {

  plog_impl_->plogger_->error(message, args);

}

template<typename M, typename... Args> void kgl::Logger::critical(M& message, Args... args) noexcept {

  plog_impl_->plogger_->critical(message, args);

}


int test_logging();


void async_example();
void syslog_example();
void android_example();
void user_defined_example();
void err_handler_example();




int test_logging()
{
  try
  {
    // Console logger with color
    auto console = spd::stdout_color_mt("console");
    console->info("Welcome to spdlog!");
    console->error("Some error message with arg{}..", 1);

    // Conditional logging example
    console->info_if(true, "Welcome to spdlog conditional logging!");

    // Formatting examples
    console->warn("Easy padding in numbers like {:08d}", 12);
    console->critical("Support for int: {0:d};  hex: {0:x};  oct: {0:o}; bin: {0:b}", 42);
    console->info("Support for floats {:03.2f}", 1.23456);
    console->info("Positional args are {1} {0}..", "too", "supported");
    console->info("{:<30}", "left aligned");

    SPDLOG_DEBUG_IF(console, true, "This is a debug log");


    spd::get("console")->info("loggers can be retrieved from a global registry using the spdlog::get(logger_name) function");


    // Create basic file logger (not rotated)
    auto my_logger = spd::basic_logger_mt("basic_logger", "logs/basic");
    my_logger->info("Some log message");

    // Create a file rotating logger with 5mb size max and 3 rotated files
    auto rotating_logger = spd::rotating_logger_mt("some_logger_name", "logs/mylogfile", 1048576 * 5, 3);
    for (int i = 0; i < 10; ++i)
      rotating_logger->info("{} * {} equals {:>10}", i, i, i*i);

    // Create a daily logger - a new file is created every day on 2:30am
    auto daily_logger = spd::daily_logger_mt("daily_logger", "logs/daily", 2, 30);
    // trigger flush if the log severity is error or higher
    daily_logger->flush_on(spd::level::err);
    daily_logger->info(123.44);

    // Customize msg format for all messages
    spd::set_pattern("*** [%H:%M:%S %z] [thread %t] %v ***");
    rotating_logger->info("This is another message with custom format");


    // Runtime log levels
    spd::set_level(spd::level::info); //Set global log level to info
    console->debug("This message shold not be displayed!");
    console->set_level(spd::level::debug); // Set specific logger's log level
    console->debug("This message shold be displayed..");

    // Compile time log levels
    // define SPDLOG_DEBUG_ON or SPDLOG_TRACE_ON
    SPDLOG_TRACE(console, "Enabled only #ifdef SPDLOG_TRACE_ON..{} ,{}", 1, 3.23);
    SPDLOG_DEBUG(console, "Enabled only #ifdef SPDLOG_DEBUG_ON.. {} ,{}", 1, 3.23);
    SPDLOG_DEBUG_IF(console, true, "This is a debug log");


    // Asynchronous logging is very fast..
    // Just call spdlog::set_async_mode(q_size) and all created loggers from now on will be asynchronous..
    async_example();

    // syslog example. linux/osx only
    syslog_example();

    // android example. compile with NDK
    android_example();

    // Log user-defined types example
    user_defined_example();

    // Change default log error handler
    err_handler_example();

    // Apply a function on all registered loggers
    spd::apply_all([&](std::shared_ptr<spdlog::logger> l)
                   {
                     l->info("End of example.");
                   });

    // Release and close all loggers
    spdlog::drop_all();
  }
    // Exceptions will only be thrown upon failed logger or sink construction (not during logging)
  catch (const spd::spdlog_ex& ex)
  {
    std::cout << "Log init failed: " << ex.what() << std::endl;
    return 1;
  }
}

void async_example()
{
  size_t q_size = 4096; //queue size must be power of 2
  spdlog::set_async_mode(q_size);
  auto async_file = spd::daily_logger_st("async_file_logger", "logs/async_log");

  for (int i = 0; i < 100; ++i)
    async_file->info("Async message #{}", i);
}

//syslog example (linux/osx/freebsd)
void syslog_example()
{
#ifdef SPDLOG_ENABLE_SYSLOG
  std::string ident = "spdlog-example";
    auto syslog_logger = spd::syslog_logger("syslog", ident, LOG_PID);
    syslog_logger->warn("This is warning that will end up in syslog.");
#endif
}

// Android example
void android_example()
{
#if defined(__ANDROID__)
  std::string tag = "spdlog-android";
    auto android_logger = spd::android_logger("android", tag);
    android_logger->critical("Use \"adb shell logcat\" to view this message.");
#endif
}

// user defined types logging by implementing operator<<
struct my_type
{
  int i;
  template<typename OStream>
  friend OStream& operator<<(OStream& os, const my_type &c)
  {
    return os << "[my_type i="<<c.i << "]";
  }
};

#include "spdlog/fmt/ostr.h" // must be included
void user_defined_example()
{
  spd::get("console")->info("user defined type: {}", my_type { 14 });
}

//
//custom error handler
//
void err_handler_example()
{
  //can be set globaly or per logger(logger->set_error_handler(..))
  spdlog::set_error_handler([](const std::string& msg)
                            {
                              std::cerr << "my err handler: " << msg << std::endl;
                            });
  spd::get("console")->info("some invalid message to trigger an error {}{}{}{}", 3);
}

