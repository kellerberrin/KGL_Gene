//
// Created by kellerberrin on 26/12/17.
//

#ifndef KGL_UTILITY_H
#define KGL_UTILITY_H


#include <string>


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

// Singleton. Defines (boost) file and time utility function.

class Utility {

public:

  Utility()=delete;
  ~Utility()=delete;

  static void getElapsedTime(double &Clock, double &System, double &User);
  static std::string filePath(const std::string& file_name, const std::string& path); // utility for "path/file"
  static std::string fileExtension(const std::string& file_name);
  static std::string fileName(const std::string& file_name);
  static std::string toupper(const std::string& s);
  static std::string trimWhiteSpace(const std::string& s);

private:

};



}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_UTILITY_H
