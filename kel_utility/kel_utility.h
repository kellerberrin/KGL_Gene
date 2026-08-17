// Copyright 2023 Kellerberrin
//

#ifndef KEL_UTILITY_H
#define KEL_UTILITY_H

#include "kel_string_hash.h"

#include <string>
#include <vector>
#include <optional>
#include <functional>
#include <cstdint>
#include <cctype>



namespace kellerberrin {   //  organization level namespace

// Object cannot be created, just supplies scope and visibility.
// Prefer namespace over class
namespace Utility {


  [[nodiscard]] std::string filePath(const std::string& file_name); // Input is "path/file.ext" function returns "path".
  [[nodiscard]] std::string filePath(const std::string& file_name, const std::string& path); // Utility to concatenate "path/file.ext".
  [[nodiscard]] std::string showPath(); // Utility to return a string with the current path plus file name.
  [[nodiscard]] std::string appendPath(const std::string& sub_directory, const std::string& path); // Utility to append subdirectory "path/sub_dir" (same as above).
  [[nodiscard]] bool fileExists(const std::string& file_path); // Check that a file exists at the file path
  [[nodiscard]] bool fileExistsCreate(const std::string& file_path); // Check that a file exists at the file path, creates zero sized file if not.
  [[nodiscard]] bool directoryExists(const std::string& path); // Check that the directory exists at the specified path
  [[nodiscard]] bool createDirectory(const std::string& path); // Create directory at the specified path
  [[nodiscard]] bool deleteDirectory(const std::string& path); // Recursively deletes the contents of a directory and all subdirectories.
  [[nodiscard]] bool recreateDirectory(const std::string& path); // Delete the directory and its contents and then recreate the directory.
  [[nodiscard]] bool directoryRenew(const std::string& path);  // If directory exists, recreate it, else just create it.
  [[nodiscard]] std::string fileExtension(const std::string& file_name); // Input is "path/file.ext", returns "ext"
  [[nodiscard]] std::string fileName(const std::string& file_name); // Input is path/file.ext", returns "file"
  [[nodiscard]] std::optional<std::string> getEnvironment(const std::string& env_var); // Translate a linux environment variable.
  [[nodiscard]] std::string toupper(const std::string& s); // Covert to upper case.
  [[nodiscard]] std::string trimAllWhiteSpace(const std::string &s); // Trim any whitespace in a string
  [[nodiscard]] std::string trimEndWhiteSpace(const std::string &s); // Only trim whitespace at either end of the string.
  [[nodiscard]] std::string trimLeadingWhiteSpace(const std::string &s); // Only trim whitespace at beginning of the string.
  [[nodiscard]] std::string trimAllChar(const std::string &s, char nc); // Returns a string with all nc char removed.
  [[nodiscard]] std::string findAndReplaceAll(const std::string& source, const std::string& search, const std::string& replace);
  [[nodiscard]] std::vector<std::string_view> viewTokenizer(const std::string_view& str_view, char delim); // Tokenize a string using delimiter chars, return std::string_view tokens.
  [[nodiscard]] std::vector<std::string> charTokenizer(const std::string_view& str_view, char delim); // Tokenize a string using delimiter chars, return std::string tokens.
  // Split string on encountering a specific char. Default version splits on first whitespace,
  [[nodiscard]] std::pair<std::string, std::string> firstSplitChar(const std::string_view& source, bool(* char_delim_fn)(char c) = [](char c)->bool { return std::isspace(static_cast<unsigned char>(c)) != 0; });
  [[nodiscard]] std::pair<double, double> process_mem_usage(); // pair.first is process vm_usage, pair.second is resident memory set.
  [[nodiscard]] std::pair<double, double> process_time_usage(); // pair.first is system CPU time usage (seconds), pair.second is user CPU time usage (seconds).
  [[nodiscard]] std::pair<double, double> process_mem_usage2(); // pair.first is process vm_usage, pair.second is physical memory used.
  [[nodiscard]] std::pair<double, double> stddev(const std::vector<double> &vec); // mean = pair.first,  sample stdev = pair.second

  // Constexpr hash algorithm for a string, can be used to hash constant strings in 'switch' statements.
  [[nodiscard]] constexpr uint64_t hash64(const std::string_view& sv) { return StringHash64::FNV_1A(sv); }


};


}   // namespace kellerberrin

#endif //KEL_UTILITY_H
