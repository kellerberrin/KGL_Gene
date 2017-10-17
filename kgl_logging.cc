//
// Created by kellerberrin on 6/09/17.
//


#include <iostream>

#include "kgl_logging.h"


namespace kgl = kellerberrin::genome;

kgl::Logger::Logger(const std::string& module, const std::string& log_file) {

  SetFormat(SPDLOG_DEFAULT_FORMAT);
  SetLevel(Severity::Trace);
  std::vector<spdlog::sink_ptr> sinks;
  sinks.push_back(std::make_shared<spdlog::sinks::stderr_sink_mt>());
  sinks.push_back(std::make_shared<spdlog::sinks::simple_file_sink_mt>(log_file));
  plog_impl_ = std::make_unique<spdlog::logger>(module, sinks.begin(), sinks.end());

}

void kgl::Logger::SetLevel(Severity level) noexcept {

  switch(level) {

    case Severity::Trace:
      spdlog::set_level(spdlog::level::trace);
      break;

    case Severity::Info:
      spdlog::set_level(spdlog::level::info);
      break;

    case Severity::Warn:
      spdlog::set_level(spdlog::level::warn);
      break;

    case Severity::Error:
      spdlog::set_level(spdlog::level::err);
      break;

    case Severity::Critical:
      spdlog::set_level(spdlog::level::critical);
      break;

  }

}

void kgl::Logger::SetFormat(const std::string& log_format) noexcept {

  spdlog::set_pattern(log_format);

}


