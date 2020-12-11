//
// Created by kellerberrin on 26/12/17.
//

#ifndef KGL_UTILITY_H
#define KGL_UTILITY_H

#include <string>
#include <vector>


namespace kellerberrin {   //  organization level namespace

// Singleton. Defines (boost) file and time utility functions and misc string functions.

class Utility {

public:

  Utility()=delete;
  ~Utility()=delete;

  [[nodiscard]] static std::string filePath(const std::string& file_name); // returns "path""
  [[nodiscard]] static std::string filePath(const std::string& file_name, const std::string& path); // utility for "path/file"
  [[nodiscard]] static bool fileExists(const std::string& file_path); // Check that a file exists at the file path
  [[nodiscard]] static std::string fileExtension(const std::string& file_name);
  [[nodiscard]] static std::string fileName(const std::string& file_name);
  [[nodiscard]] static std::string toupper(const std::string& s);
  [[nodiscard]] static std::string trimAllWhiteSpace(const std::string &s); // Trim any whitespace in a string
  [[nodiscard]] static std::string trimEndWhiteSpace(const std::string &s); // Only trim whitespace at either end of the string.
  [[nodiscard]] static std::string trimAllChar(const std::string &s, const char nc); // Returns a string with all nc char removed.
  [[nodiscard]] static std::string findAndReplaceAll(const std::string& source, const std::string& search, const std::string& replace);
  [[nodiscard]] static std::vector<std::string> tokenizer(const std::string& str, const std::string& delims = ",");
  [[nodiscard]] static std::vector<std::string_view> view_tokenizer(const std::string_view& str_view, const char delim);
  [[nodiscard]] static std::vector<std::string> char_tokenizer(const std::string& str, const char delim);
  [[nodiscard]] static std::pair<double, double> process_mem_usage();
  [[nodiscard]] static std::pair<double, double> stddev(const std::vector<double> &vec); // mean = first,  stdev = second
  static void getElapsedTime(double &Clock, double &System, double &User);

private:

};



}   // end namespace

#endif //KGL_UTILITY_H
