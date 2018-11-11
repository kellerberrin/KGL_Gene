//
// Created by kellerberrin on 26/12/17.
//

#ifndef KGL_UTILITY_H
#define KGL_UTILITY_H


#include <string>


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

// Singleton. Defines (boost) file and time utility functions and misc string functions.

class Utility {

public:

  Utility()=delete;
  ~Utility()=delete;

  static void getElapsedTime(double &Clock, double &System, double &User);
  static std::string filePath(const std::string& file_name, const std::string& path); // utility for "path/file"
  static std::string fileExtension(const std::string& file_name);
  static std::string fileName(const std::string& file_name);
  static std::string toupper(const std::string& s);
  static std::string trimAllWhiteSpace(const std::string &s); // Trim any whitespace in a string
  static std::string trimEndWhiteSpace(const std::string &s); // Only trim whitespace at either end of the string.

private:

};



}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_UTILITY_H
