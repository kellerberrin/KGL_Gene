//
// Created by kellerberrin on 26/12/17.
//


#include <iostream>
#include "kgl_utility.h"
#include <boost/timer/timer.hpp>
#define BOOST_FILESYSTEM_NO_DEPRECATED // Recommended by boost filesystem documentation.
#include <boost/filesystem.hpp>


// Define namespace alias
namespace bt = boost::timer;
namespace fs = boost::filesystem;
namespace kgl = kellerberrin::genome;



// Hide the boost cpu timer in an anonymous namespace.
namespace {  bt::cpu_timer cpu_timer; }

void kgl::Utility::getElapsedTime(double &Clock, double &User, double &System) {

  Clock = 0; User = 0; System = 0;
  bt::cpu_times elapsedtime = cpu_timer.elapsed();
  Clock = elapsedtime.wall / 1e09; // Convert from nanoseconds to seconds
  User = elapsedtime.user / 1e09;
  System = elapsedtime.system / 1e09;

}


// Returns the filename with the path directory appended to it "path/file".
std::string kgl::Utility::filePath(const std::string& file_name, const std::string& path) {

  fs::path directory_path = fs::path(path);
  fs::path file_path = directory_path / fs::path(file_name);
  return file_path.string();

}

// Returns the filename extension.
std::string kgl::Utility::fileExtension(const std::string& file_name) {

  fs::path file_path(file_name);
  return file_path.extension().string();  // returns the extension with a '.' e.g. "example.txt" returns ".txt"

}

// Returns the filename extension.
std::string kgl::Utility::fileName(const std::string& file_name) {

  fs::path file_path(file_name);
  return file_path.filename().string();  // returns the file name without path. "/path/file.ext" returns "file.ext".

}

// Returns uppercase string.
std::string kgl::Utility::toupper(const std::string& s) {

  std::string upper_string;
  auto lambda_to_upper = [](unsigned char c){ return std::toupper(c); };
  std::transform(s.begin(), s.end(), std::back_inserter(upper_string), lambda_to_upper);

  return upper_string;

}

// Returns trimmed string.
std::string kgl::Utility::trimWhiteSpace(const std::string& s) {

  std::string clean_string;
  auto lambda_not_whitespace = [](unsigned char c){ return not std::isspace(c); };
  std::copy_if(s.begin(), s.end(), back_inserter(clean_string), lambda_not_whitespace);

  return clean_string;

}
