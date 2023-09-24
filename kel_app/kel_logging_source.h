//
// Created by kellerberrin on 24/09/23.
//

#ifndef KEL_LOGGING_SOURCE_H
#define KEL_LOGGING_SOURCE_H


#include <spdlog/spdlog.h>  // Implement the logger using the spdlog library

#include <string_view>
#include <memory>
#include <iostream>
#include <source_location>


namespace kellerberrin {   //  organization level namespace


[[nodiscard]] constexpr auto getSourceLocation(const std::source_location &location) {

  return spdlog::source_loc{location.file_name(),
                            static_cast<std::int32_t>(location.line()),
                            location.function_name()};

}

struct FormatLocation {

  std::string_view format_view;
  spdlog::source_loc spdlog_locaton;

  template<typename String>
  FormatLocation(const String &format_string,
                 const std::source_location &location = std::source_location::current())
      : format_view{format_string}, spdlog_locaton{getSourceLocation(location)} {}

};

template<typename... Args>
void warn(FormatLocation format_location, Args &&...args) {

  spdlog::default_logger_raw()->log(format_location.spdlog_locaton,
                                    spdlog::level::warn,
                                    format_location.format_view,
                                    std::forward<Args>(args)...);

}

/*
int main() {
  spdlog::set_level(spdlog::level::trace);
  std::string message = "hello {}\n";
  auto value = 42;
  logging::warn(message, value);
}
*/


} // Namespace

#endif //KEL_LOGGING_SOURCE_H
