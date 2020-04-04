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

  static void getElapsedTime(double &Clock, double &System, double &User);
  [[nodiscard]] static std::string filePath(const std::string& file_name, const std::string& path); // utility for "path/file"
  [[nodiscard]] static bool fileExists(const std::string& file_path); // Check that a file exists at the file path
  [[nodiscard]] static std::string fileExtension(const std::string& file_name);
  [[nodiscard]] static std::string fileName(const std::string& file_name);
  [[nodiscard]] static std::string toupper(const std::string& s);
  [[nodiscard]] static std::string trimAllWhiteSpace(const std::string &s); // Trim any whitespace in a string
  [[nodiscard]] static std::string trimEndWhiteSpace(const std::string &s); // Only trim whitespace at either end of the string.
  [[nodiscard]] static std::string findAndReplaceAll(const std::string& source, const std::string& search, const std::string& replace);
  [[nodiscard]] static std::vector<std::string> tokenizer(const std::string& str, const std::string& delims = ",");


private:

};



}   // end namespace

#endif //KGL_UTILITY_H
