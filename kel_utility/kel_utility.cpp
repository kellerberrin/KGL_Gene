//
// Created by kellerberrin on 26/12/17.
//


#include "kel_utility.h"
#include "kel_exec_env.h"

#include <boost/timer/timer.hpp>

#include <unistd.h>
#include <ios>
#include <string>
#include <numeric>
#include <fstream>
#include <filesystem>


namespace bt = boost::timer;
namespace kel = kellerberrin;
namespace fs = std::filesystem;


// Given a file name as "path/file", returns "path"
std::string kel::Utility::filePath(const std::string& file_name) {

  fs::path file_path = fs::path(file_name);
  fs::path directory = file_path.parent_path();
  return directory.string();

}


// Returns the filename with the path directory appended to it "path/file".
std::string kel::Utility::filePath(const std::string& file_name, const std::string& path) {

  fs::path directory_path = fs::path(path);
  fs::path file_path = directory_path / fs::path(file_name);
  return file_path.string();

}

// Check that a file exists at the file path.
bool kel::Utility::fileExists(const std::string& file_path) {

  bool file_exists = fs::exists(fs::path(file_path));
  return file_exists;

}

// Check that a file exists at the file path, creates zero sized file if not.
bool kel::Utility::fileExistsCreate(const std::string& file_path) {

  if (not fileExists(file_path)) {

    std::ofstream empty_file(file_path);
    return empty_file.good();

  }

  return true;

}


// Check that the directory exists at the specified path
bool kel::Utility::directoryExists(const std::string& path) {

  return fileExists(path);

}

// Create directory at the specified path, returns true if directory already exists.
bool kel::Utility::createDirectory(const std::string& path) {

  if (directoryExists(path)) {

    return true;

  }
  return fs::create_directory(fs::path(path));

}

// Recursively delete the contents of a directory and then the directory.
bool kel::Utility::deleteDirectory(const std::string& path) {

  return fs::remove_all(fs::path(path));

}

// Delete the directory and its contents and then recreate the directory.
bool kel::Utility::recreateDirectory(const std::string& path) {

  bool result = deleteDirectory(path);
  if (result) {

    result = createDirectory(path);

  }

  return result;

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


std::optional<std::string> kel::Utility::getEnvironment(const std::string& env_var) {

  const char *val = std::getenv(env_var.c_str());

  if (val == nullptr ) {

    return std::nullopt;

  }
  else {

    return std::string(val);

  }

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

// Returns a string with all nc char removed.
std::string kel::Utility::trimAllChar(const std::string &s, char nc) {

  std::string mod_string;
  auto lambda_not_char = [nc](unsigned char c){ return c != nc; };
  std::copy_if(s.begin(), s.end(), back_inserter(mod_string), lambda_not_char);

  return mod_string;

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


std::vector<std::string_view> kel::Utility::viewTokenizer(const std::string_view& str_view, char delim) {

  std::vector<std::string_view> token_vector;
  size_t token_index{0};
  size_t index{0};

  for (; index < str_view.size(); ++index) {

    if (str_view[index] == delim) {

      token_vector.emplace_back(&str_view[token_index], (index-token_index));
      token_index = index + 1;

    }

  }

  if (token_index > index) {

    token_vector.emplace_back();

  } else {

    token_vector.emplace_back(&str_view[token_index], (index-token_index));

  }

  return token_vector;

}


std::vector<std::string> kel::Utility::charTokenizer(const std::string_view& str_view, char delim) {

  std::vector<std::string_view> view_vector = viewTokenizer(str_view, delim);
  std::vector<std::string> str_vector;
  str_vector.reserve(view_vector.size());
  for (auto const& view : view_vector) {

    str_vector.emplace_back(view);

  }

  return str_vector;

}


// Split string on encountering char. Default version splits on first whitespace
std::pair<std::string, std::string> kel::Utility::firstSplit(const std::string& source, bool(* char_delim_fn)(char c)) {

  std::string first;
  std::string second;
  bool split{false};

  for (const char c : source) {

    if (char_delim_fn(c)) {

      split = true;

    }

    if (not split) {

      first += c;

    } else {

      second += c;

    }

  }

  return {first, second};

}


//////////////////////////////////////////////////////////////////////////////
//

// Attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0
std::pair<double, double> kel::Utility::process_mem_usage()
{

  // 'file' stat seems to give the most reliable results
  static const char* STAT_STREAM = "/proc/self/stat";


  std::ifstream stat_stream(STAT_STREAM,std::ios_base::in);

  if (not stat_stream.good()) {

    ExecEnv::log().error("Utility::process_mem_usage; Unable to open file stats: {}", STAT_STREAM);
    return {0.0, 0.0};

  }

  // Dummy vars for leading entries in stat that we don't care about.
  //
  std::string pid, comm, state, ppid, pgrp, session, tty_nr;
  std::string tpgid, flags, minflt, cminflt, majflt, cmajflt;
  std::string utime, stime, cutime, cstime, priority, nice;
  std::string _o_, itrealvalue, starttime;

  // The two fields we want
  //
  unsigned long vsize;
  long rss;

  stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
              >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
              >> utime >> stime >> cutime >> cstime >> priority >> nice
              >> _o_ >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

  stat_stream.close();

  // Calc to KBytes.
  long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
  double vm_usage     = vsize / 1024.0;
  double resident_set = rss * page_size_kb;

  // Convert to GBytes.
  vm_usage = vm_usage / (1024.0 * 1024.0);
  resident_set = resident_set / (1024.0 * 1024.0);

  return std::pair<double, double>(vm_usage, resident_set);

}



//////////////////////////////////////////////////////////////////////////////
//
// pair.first is process vm_usage, pair.second is physical memory used.
// On failure, returns 0.0, 0.0
std::pair<double, double> kel::Utility::process_mem_usage2()
{

  // The status pseudo-file.
  static const char* STAT_STREAM = "/proc/self/status";
  // vm used.
  static const char* VM_USAGE = "VmSize:";
  // Physical mem used.
  static const char* PM_USAGE = "VmRSS:";
  // Trimmed text.
  static const char* KB_TRIM = "kB";
  // Scale to MegaBytes.
  static const double SCALE_MB = 1024.0;

  std::ifstream stat_stream(STAT_STREAM,std::ios_base::in);

  if (not stat_stream.good()) {

    ExecEnv::log().error("Utility::process_mem_usage2; Unable to open file stats: {}", STAT_STREAM);
    return {0.0, 0.0};

  }

  std::pair<double, double> mem_stats{0.0, 0.0};
  std::string stat_line;
  while (std::getline(stat_stream, stat_line)) {

    if (stat_line.find(VM_USAGE) != std::string::npos) {

      std::string trimmed_stat = Utility::findAndReplaceAll(stat_line, VM_USAGE, "");
      trimmed_stat = Utility::findAndReplaceAll(trimmed_stat, KB_TRIM, "");
      trimmed_stat = Utility::trimAllWhiteSpace(trimmed_stat);
      try {

        mem_stats.first = std::stod(trimmed_stat);

      } catch(...) {

        ExecEnv::log().error("process_mem_usage2; cannot process {} with value: {}", VM_USAGE, stat_line);

      }

    }

    if (stat_line.find(PM_USAGE) != std::string::npos) {

      std::string trimmed_stat = Utility::findAndReplaceAll(stat_line, PM_USAGE, "");
      trimmed_stat = Utility::findAndReplaceAll(trimmed_stat, KB_TRIM, "");
      trimmed_stat = Utility::trimAllWhiteSpace(trimmed_stat);
      try {

        mem_stats.second = std::stod(trimmed_stat);

      } catch(...) {

        ExecEnv::log().error("process_mem_usage2; cannot process {} with value: {}", PM_USAGE, stat_line);

      }

    }

  }

  return { mem_stats.first / SCALE_MB, mem_stats.second / SCALE_MB };

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
  auto variance_func = [&mean](const double& accumulator, const double& val)->double
  {
    return accumulator + ((val - mean) * (val - mean));
  };

  // Sample variance
  double variance = std::accumulate(vec.begin(), vec.end(), 0.0, variance_func) / (size - 1.0);

  return {mean, std::sqrt(variance) };
}


// Hide the boost cpu timer in an anonymous namespace.
namespace {  bt::cpu_timer cpu_timer; }

void kel::Utility::getRuntime(double &Clock, double &User, double &System) {

  Clock = 0; User = 0; System = 0;
  bt::cpu_times elapsedtime = cpu_timer.elapsed();
  Clock = elapsedtime.wall / 1e09; // Convert from nanoseconds to seconds
  User = elapsedtime.user / 1e09;
  System = elapsedtime.system / 1e09;

}


