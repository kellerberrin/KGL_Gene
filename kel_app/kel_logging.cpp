// Copyright 2023 Kellerberrin
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
// documentation files (the "Software"), to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
// and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
// WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
// OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//
//


#include <iostream>

#include "kel_logging.h"

#include "spdlog/sinks/basic_file_sink.h"


namespace kel = kellerberrin;




kel::Logger::Logger(const std::string& module, const std::string& log_file) {

  setFormat(SPDLOG_DEFAULT_FORMAT);
  setLevel(Severity::Trace);
  std::vector<spdlog::sink_ptr> sinks;
  sinks.push_back(std::make_shared<spdlog::sinks::ansicolor_stderr_sink_mt>());
  sinks.push_back(std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file));
  log_impl_ptr_ = std::make_unique<spdlog::logger>(module, sinks.begin(), sinks.end());

}


void kel::Logger::setLevel(Severity level) noexcept {

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

void kel::Logger::setFormat(const std::string& log_format) noexcept {

  spdlog::set_pattern(log_format);

}



