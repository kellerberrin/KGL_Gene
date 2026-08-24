// Copyright 2023 Kellerberrin
//

#ifndef KEL_UTILITY_H
#define KEL_UTILITY_H

#include "kel_string_hash.h"

#include <string>
#include <string_view>
#include <vector>
#include <optional>
#include <utility>
#include <functional>
#include <cstdint>
#include <cctype>



namespace kellerberrin {   //  organization level namespace

// Object cannot be created, just supplies scope and visibility.
// Prefer namespace over class
namespace Utility {


  /// Returns the directory portion of "path/file.ext" as "path".
  [[nodiscard]] std::string filePath(const std::string& file_name);
  /// Concatenates file_name onto path, returning "path/file.ext".
  [[nodiscard]] std::string filePath(const std::string& file_name, const std::string& path);
  /// Returns the current working directory.
  [[nodiscard]] std::string showPath();
  /// Appends sub_directory onto path, returning "path/sub_dir".
  [[nodiscard]] std::string appendPath(const std::string& sub_directory, const std::string& path);
  /// Returns true if a regular file exists at file_path.
  [[nodiscard]] bool fileExists(const std::string& file_path);
  /// Returns true if a file exists at file_path; creates a zero-size file if not.
  [[nodiscard]] bool fileExistsCreate(const std::string& file_path);
  /// Returns true if a directory exists at path.
  [[nodiscard]] bool directoryExists(const std::string& path);
  /// Creates a directory at path. Returns true if it already exists or was created.
  [[nodiscard]] bool createDirectory(const std::string& path);
  /// Recursively deletes a directory and all its contents. Returns true if anything was removed.
  [[nodiscard]] bool deleteDirectory(const std::string& path);
  /// Deletes a directory and its contents, then recreates the directory.
  [[nodiscard]] bool recreateDirectory(const std::string& path);
  /// If directory exists, recreates it; otherwise creates it.
  [[nodiscard]] bool directoryRenew(const std::string& path);
  /// Returns the extension of "path/file.ext" as ".ext" (includes the dot).
  [[nodiscard]] std::string fileExtension(const std::string& file_name);
  /// Returns the filename portion of "path/file.ext" as "file.ext".
  [[nodiscard]] std::string fileName(const std::string& file_name);
  /// Translates a Linux environment variable. Returns std::nullopt if not set.
  [[nodiscard]] std::optional<std::string> getEnvironment(const std::string& env_var);
  /// Converts the string to upper case.
  [[nodiscard]] std::string toupper(const std::string& s);
  /// Removes all whitespace characters from the string.
  [[nodiscard]] std::string trimAllWhiteSpace(const std::string &s);
  /// Trims whitespace at both ends of the string.
  [[nodiscard]] std::string trimEndWhiteSpace(const std::string &s);
  /// Trims whitespace at the beginning of the string only.
  [[nodiscard]] std::string trimLeadingWhiteSpace(const std::string &s);
  /// Returns a string with all occurrences of character nc removed.
  [[nodiscard]] std::string trimAllChar(const std::string &s, char nc);
  /// Replaces all occurrences of search with replace in source.
  [[nodiscard]] std::string findAndReplaceAll(const std::string& source, const std::string& search, const std::string& replace);
  /// Tokenizes a string_view using a delimiter character, returning string_view tokens.
  [[nodiscard]] std::vector<std::string_view> viewTokenizer(const std::string_view& str_view, char delim);
  /// Tokenizes a string_view using a delimiter character, returning string tokens.
  [[nodiscard]] std::vector<std::string> charTokenizer(const std::string_view& str_view, char delim);
  /// Splits a string at the first character matching char_delim_fn. Default splits on first whitespace.
  [[nodiscard]] std::pair<std::string, std::string> firstSplitChar(const std::string_view& source, bool(* char_delim_fn)(char c) = [](char c)->bool { return std::isspace(static_cast<unsigned char>(c)) != 0; });
  /// Returns {virtual_memory_GB, resident_set_GB} from /proc/self/stat. Returns {0.0, 0.0} on failure.
  [[nodiscard]] std::pair<double, double> process_mem_usage();
  /// Returns {system_cpu_seconds, user_cpu_seconds} from /proc/self/stat. Returns {0.0, 0.0} on failure.
  [[nodiscard]] std::pair<double, double> process_time_usage();
  /// Returns {virtual_memory_MB, physical_memory_MB} from /proc/self/status. Returns {0.0, 0.0} on failure.
  [[nodiscard]] std::pair<double, double> process_mem_usage2();
  /// Returns {mean, sample_standard_deviation} of the input vector.
  [[nodiscard]] std::pair<double, double> stddev(const std::vector<double> &vec);

  /// Constexpr hash algorithm for a string, can be used to hash constant strings in 'switch' statements.
  [[nodiscard]] constexpr uint64_t hash64(const std::string_view& sv) { return StringHash64::FNV_1A(sv); }


};


}   // namespace kellerberrin

#endif //KEL_UTILITY_H
