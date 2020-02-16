//
// Created by kellerberrin on 26/12/17.
//


#include <iostream>
#include "kel_utility.h"
#include <boost/timer/timer.hpp>
#define BOOST_FILESYSTEM_NO_DEPRECATED // Recommended by boost filesystem documentation.
#include <boost/filesystem.hpp>


// Define namespace alias
namespace bt = boost::timer;
namespace fs = boost::filesystem;
namespace kel = kellerberrin;



// Hide the boost cpu timer in an anonymous namespace.
namespace {  bt::cpu_timer cpu_timer; }

void kel::Utility::getElapsedTime(double &Clock, double &User, double &System) {

  Clock = 0; User = 0; System = 0;
  bt::cpu_times elapsedtime = cpu_timer.elapsed();
  Clock = elapsedtime.wall / 1e09; // Convert from nanoseconds to seconds
  User = elapsedtime.user / 1e09;
  System = elapsedtime.system / 1e09;

}


// Returns the filename with the path directory appended to it "path/file".
std::string kel::Utility::filePath(const std::string& file_name, const std::string& path) {

  fs::path directory_path = fs::path(path);
  fs::path file_path = directory_path / fs::path(file_name);
  return file_path.string();

}

// Check that a file exists at the file path.
bool kel::Utility::fileExists(const std::string& file_path) {

  boost::system::error_code error_code;
  bool file_exists = fs::exists(fs::path(file_path), error_code);

  if (error_code.value() != boost::system::errc::success) {

    return false;

  }

  return file_exists;

}

// Returns the filename extension.
std::string kel::Utility::fileExtension(const std::string& file_name) {

  fs::path file_path(file_name);
  return file_path.extension().string();  // returns the extension with a '.' e.g. "example.txt" returns ".txt"

}

// Returns the filename extension.
std::string kel::Utility::fileName(const std::string& file_name) {

  fs::path file_path(file_name);
  return file_path.filename().string();  // returns the file name without path. "/path/file.ext" returns "file.ext".

}

// Returns uppercase string.
std::string kel::Utility::toupper(const std::string& s) {

  std::string upper_string;
  auto lambda_to_upper = [](unsigned char c){ return std::toupper(c); };
  std::transform(s.begin(), s.end(), std::back_inserter(upper_string), lambda_to_upper);

  return upper_string;

}

// Returns trimmed string.
std::string kel::Utility::trimAllWhiteSpace(const std::string &s) {

  std::string clean_string;
  auto lambda_not_whitespace = [](unsigned char c){ return not std::isspace(c); };
  std::copy_if(s.begin(), s.end(), back_inserter(clean_string), lambda_not_whitespace);

  return clean_string;

}

// Only trim whitespace at either end of the string.
std::string kel::Utility::trimEndWhiteSpace(const std::string &s) {

  std::string start_trimmed_string;

  auto it = s.begin();
  while (it != s.end()) {

    if (not std::isspace(*it)) {

      break;

    }

    ++it;

  }

  while(it != s.end()) {

    start_trimmed_string.push_back(*it);
    ++it;

  }

  std::string trimmed_string;

  auto rit = start_trimmed_string.rbegin();

  while(rit != start_trimmed_string.rend()) {

    if (not std::isspace(*rit)) {

      break;

    }

    ++rit;

  }

  while(rit != start_trimmed_string.rend()) {

    trimmed_string.push_back(*rit);
    ++rit;

  }

  std::reverse(trimmed_string.begin(), trimmed_string.end());

  return trimmed_string;

}


std::string kel::Utility::findAndReplaceAll(const std::string& source, const std::string& search, const std::string& replace)
{

  std::string modified = source;
  // Get the first occurrence
  size_t pos = modified.find(search);

  // Repeat till end is reached
  while( pos != std::string::npos)
  {
    // Replace this occurrence of Sub String
    modified.replace(pos, search.size(), replace);
    // Get the next occurrence from the position after the replaced string.
    pos = modified.find(search, pos + replace.size());

  }

  return modified;

}