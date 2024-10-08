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


namespace kellerberrin {   //  organization level namespace

// Object cannot be created, just supplies scope and visibility.
// Defines file and time utility functions and misc string/parsing functions.

class Utility {

public:

  Utility()=delete;
  ~Utility()=delete;

  [[nodiscard]] static std::string filePath(const std::string& file_name); // Input is "path/file.ext" function returns "path".
  [[nodiscard]] static std::string filePath(const std::string& file_name, const std::string& path); // Utility to concatenate "path/file.ext".
  [[nodiscard]] static std::string appendPath(const std::string& sub_directory, const std::string& path); // Utility to append subdirectory "path/sub_dir" (same as above).
  [[nodiscard]] static bool fileExists(const std::string& file_path); // Check that a file exists at the file path
  [[nodiscard]] static bool fileExistsCreate(const std::string& file_path); // Check that a file exists at the file path, creates zero sized file if not.
  [[nodiscard]] static bool directoryExists(const std::string& path); // Check that the directory exists at the specified path
  [[nodiscard]] static bool createDirectory(const std::string& path); // Create directory at the specified path
  [[nodiscard]] static bool deleteDirectory(const std::string& path); // Recursively deletes the contents of a directory and all subdirectories.
  [[nodiscard]] static bool recreateDirectory(const std::string& path); // Delete the directory and its contents and then recreate the directory.
  [[nodiscard]] static bool directoryRenew(const std::string& path);  // If directory exists, recreate it, else just create it.
  [[nodiscard]] static std::string fileExtension(const std::string& file_name); // Input is "path/file.ext", returns "ext"
  [[nodiscard]] static std::string fileName(const std::string& file_name); // Input is path/file.ext", returns "file"
  [[nodiscard]] static std::optional<std::string> getEnvironment(const std::string& env_var); // Translate a linux environment variable.
  [[nodiscard]] static std::string toupper(const std::string& s); // Covert to upper case.
  [[nodiscard]] static std::string trimAllWhiteSpace(const std::string &s); // Trim any whitespace in a string
  [[nodiscard]] static std::string trimEndWhiteSpace(const std::string &s); // Only trim whitespace at either end of the string.
  [[nodiscard]] static std::string trimAllChar(const std::string &s, char nc); // Returns a string with all nc char removed.
  [[nodiscard]] static std::string findAndReplaceAll(const std::string& source, const std::string& search, const std::string& replace);
  [[nodiscard]] static std::vector<std::string_view> viewTokenizer(const std::string_view& str_view, char delim); // Tokenize a string using delimiter chars, return std::string_view tokens.
  [[nodiscard]] static std::vector<std::string> charTokenizer(const std::string_view& str_view, char delim); // Tokenize a string using delimiter chars, return std::string tokens.
  // Split string on encountering char. Default version splits on first whitespace,
  [[nodiscard]] static std::pair<std::string, std::string> firstSplit(const std::string& source, bool(* char_delim_fn)(char c) = [](char c)->bool { return std::isspace(c) != 0; });
  [[nodiscard]] static std::pair<double, double> process_mem_usage(); // pair.first is process vm_usage, pair.second is resident memory set.
  [[nodiscard]] static std::pair<double, double> process_time_usage(); // pair.first is system CPU time usage (seconds), pair.second is user CPU time usage (seconds).
  [[nodiscard]] static std::pair<double, double> process_mem_usage2(); // pair.first is process vm_usage, pair.second is physical memory used.
  [[nodiscard]] static std::pair<double, double> stddev(const std::vector<double> &vec); // mean = pair.first,  sample stdev = pair.second
  // Constexpr hash algorithm for a string, can be used to hash constant strings in 'switch' statements.
  [[nodiscard]] constexpr static uint32_t hash(const std::string_view& sv) { return UtilityStringHash::crcStringHash(sv); }

private:

};



}   // end namespace

#endif //KEL_UTILITY_H
