//
// Created by kellerberrin on 26/12/17.
//


#include "kel_utility.h"
#include "kel_exec_env.h"

#define BOOST_FILESYSTEM_NO_DEPRECATED 1 // Recommended by boost filesystem documentation.
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <boost/timer/timer.hpp>

#include <algorithm>
#include <unistd.h>
#include <ios>
#include <string>
#include <numeric>


namespace fs = boost::filesystem;
namespace bt = boost::timer;
namespace kel = kellerberrin;



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


std::vector<std::string> kel::Utility::tokenizer(const std::string& str, const std::string& delims) {

  if (delims.size() == 1) {

    return char_tokenizer(str, delims[0]);

  }

  std::vector<std::string> token_vector;
  boost::tokenizer<boost::char_separator<char>> tokens(str, boost::char_separator<char>(delims.c_str()));
  std::copy(tokens.begin(), tokens.end(), std::back_inserter(token_vector));

  return token_vector;

}


std::vector<std::string_view> kel::Utility::view_tokenizer(const std::string_view& str_view, const char delim) {

  std::vector<std::string_view> token_vector;
  size_t token_index{0};
  size_t token_count{0};

  for (size_t index = 0; index < str_view.size(); ++index) {

    if (str_view[index] == delim) {

      token_vector.emplace_back(str_view.substr(token_index, token_count));
      token_index = index + 1;
      token_count = 0;

    } else {

      ++token_count;

    }

  }

  if (token_index + token_count != str_view.size()) {

    ExecEnv::log().error("Utility::view_tokenizer, final token index: {}, token count: {} does not equal string_view size: {}",
                         token_index, token_count, str_view.size());

  } else {

    if (token_count == 0) {

      token_vector.emplace_back(std::string_view());

    } else {

      token_vector.emplace_back(str_view.substr(token_index, token_count));

    }

  }

  return token_vector;

}


std::vector<std::string> kel::Utility::char_tokenizer(const std::string& str, const char delim) {

  std::vector<std::string_view> view_vector = view_tokenizer(str, delim);
  std::vector<std::string> str_vector;
  str_vector.reserve(view_vector.size());
  for (auto const& view : view_vector) {

    str_vector.emplace_back(std::string(view));

  }

  return str_vector;

}

//////////////////////////////////////////////////////////////////////////////
//
// pair.first is process vm_usage, pair.second is resident memory set.
// Attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0
std::pair<double, double> kel::Utility::process_mem_usage()
{

  using std::ios_base;
  using std::ifstream;
  using std::string;
  double vm_usage{0.0};
  double resident_set{0.0};

  // 'file' stat seems to give the most reliable results
  //
  ifstream stat_stream("/proc/self/stat",ios_base::in);

  // dummy vars for leading entries in stat that we don't care about
  //
  string pid, comm, state, ppid, pgrp, session, tty_nr;
  string tpgid, flags, minflt, cminflt, majflt, cmajflt;
  string utime, stime, cutime, cstime, priority, nice;
  string O, itrealvalue, starttime;

  // the two fields we want
  //
  unsigned long vsize;
  long rss;

  stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
              >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
              >> utime >> stime >> cutime >> cstime >> priority >> nice
              >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

  stat_stream.close();

  // Calc to KBytes.
  long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
  vm_usage     = vsize / 1024.0;
  resident_set = rss * page_size_kb;

  // Convert to GBytes.
  vm_usage = vm_usage / (1024.0 * 1024.0);
  resident_set = resident_set / (1024.0 * 1024.0);

  return std::pair<double, double>(vm_usage, resident_set);

}

// pair first is the mean, second is the sample standard deviation
std::pair<double, double> kel::Utility::stddev(const std::vector<double> &vec)
{
  if (vec.empty()) {

    return {0.0, 0.0};

  }

  size_t sz = vec.size();
  if (sz == 1) {

    return {vec.front(), 0.0};

  }

  double size = static_cast<double>(sz);

  // Calculate the mean
   double mean = std::accumulate(vec.begin(), vec.end(), 0.0) / size;

  // Now calculate the variance
  auto variance_func = [&mean](double accumulator, const double& val)
  {
    return accumulator + ((val - mean) * (val - mean));
  };

  // Sample variance
  double variance = std::accumulate(vec.begin(), vec.end(), 0.0, variance_func) / (size - 1.0);

  return {mean, std::sqrt(variance) };
}


// Hide the boost cpu timer in an anonymous namespace.
namespace {  bt::cpu_timer cpu_timer; }

void kel::Utility::getElapsedTime(double &Clock, double &User, double &System) {

  Clock = 0; User = 0; System = 0;
  bt::cpu_times elapsedtime = cpu_timer.elapsed();
  Clock = elapsedtime.wall / 1e09; // Convert from nanoseconds to seconds
  User = elapsedtime.user / 1e09;
  System = elapsedtime.system / 1e09;

}

