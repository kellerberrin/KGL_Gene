// Copyright 2023 Kellerberrin
//


#include "kel_utility.h"
#include "kel_exec_env.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <unistd.h>
#include <vector>

namespace kel = kellerberrin;
namespace fs = std::filesystem;

#define  UTILITY_USE_NEW 1
#ifdef UTILITY_USE_NEW

// ---------------------------------------------------------------------------
// Path / filesystem helpers
// ---------------------------------------------------------------------------

std::string kel::Utility::filePath(const std::string& file_name) {

  return fs::path(file_name).parent_path().string();

}

std::string kel::Utility::filePath(const std::string& file_name, const std::string& path) {

  return (fs::path(path) / file_name).string();

}

// NOTE: Consider renaming this to currentPath() — "showPath" does not describe
//       that it returns the current working directory.
std::string kel::Utility::showPath() {

  return fs::current_path().string();

}

std::string kel::Utility::appendPath(const std::string& sub_directory, const std::string& path) {

  return filePath(sub_directory, path);

}

bool kel::Utility::fileExists(const std::string& file_path) {

  std::error_code ec;
  bool exists = fs::is_regular_file(file_path, ec);

  if (ec) {

    ExecEnv::log().error("Utility::fileExists; error checking {}: {}", file_path, ec.message());
    return false;

  }
  return exists;

}

bool kel::Utility::fileExistsCreate(const std::string& file_path) {

  if (fileExists(file_path)) {

    return true;

  }
  std::ofstream empty_file(file_path);
  if (!empty_file.is_open()) {

    ExecEnv::log().error("Utility::fileExistsCreate; failed to create {}", file_path);
    return false;

  }
  return true;

}

bool kel::Utility::directoryExists(const std::string& path) {

  return fs::is_directory(fs::path(path));

}

bool kel::Utility::createDirectory(const std::string& path) {

  if (directoryExists(path)) {

    return true;

  }

  std::error_code ec;
  bool created = fs::create_directory(path, ec);
  if (ec) {

    ExecEnv::log().error("Utility::createDirectory; failed to create {}: {}", path, ec.message());
    return false;

  }
  return created;

}

bool kel::Utility::deleteDirectory(const std::string& path_string) {
    std::error_code ec;
    const std::uintmax_t removed = fs::remove_all(path_string, ec);
    if (ec) {

      ExecEnv::log().error("Utility::deleteDirectory; failed to remove {}: {}", path_string, ec.message());
      return false;

    }
    // remove_all returns 0 if the path did not exist, which maps to "nothing removed".
    return removed > 0;
}

bool kel::Utility::recreateDirectory(const std::string& path) {

  if (!deleteDirectory(path)) {

    return false;
  }

  return createDirectory(path);

}

bool kel::Utility::directoryRenew(const std::string& path) {

  if (directoryExists(path)) {

    return recreateDirectory(path);

  }

  return createDirectory(path);

}

std::string kel::Utility::fileExtension(const std::string& file_name) {

  return fs::path(file_name).extension().string();

}

std::string kel::Utility::fileName(const std::string& file_name) {

  return fs::path(file_name).filename().string();

}

// ---------------------------------------------------------------------------
// Environment / string helpers
// ---------------------------------------------------------------------------

std::optional<std::string> kel::Utility::getEnvironment(const std::string& env_var) {
    if (const char* val = std::getenv(env_var.c_str())) {
        return std::string(val);
    }
    return std::nullopt;
}

std::string kel::Utility::toupper(const std::string& s) {

  std::string upper_string;
  upper_string.reserve(s.size());

  auto toupper_lambda = [](unsigned char c)->char { return static_cast<char>(std::toupper(c)); };
  std::transform(s.begin(), s.end(), std::back_inserter(upper_string), toupper_lambda);

  return upper_string;
}

std::string kel::Utility::trimAllWhiteSpace(const std::string& s) {

  std::string clean_string;
  clean_string.reserve(s.size());

  auto ws_pred = [](unsigned char c)->bool { return !std::isspace(c); };
  std::copy_if(s.begin(), s.end(), std::back_inserter(clean_string), ws_pred);
  return clean_string;

}

std::string kel::Utility::trimAllChar(const std::string& s, char nc) {

  std::string mod_string;
  mod_string.reserve(s.size());

  const auto target = static_cast<unsigned char>(nc);
  auto trim_pred = [target](unsigned char c)->bool { return c != target; };
  std::copy_if(s.begin(), s.end(), std::back_inserter(mod_string), trim_pred);

  return mod_string;

}

std::string kel::Utility::trimEndWhiteSpace(const std::string& s) {

  const auto not_space = [](unsigned char c)->bool { return !std::isspace(c); };
  const auto first = std::find_if(s.begin(), s.end(), not_space);
  if (first == s.end()) {
        return {};
  }

  const auto last = std::find_if(s.rbegin(), s.rend(), not_space).base();
  return {first, last};

}

std::string kel::Utility::trimLeadingWhiteSpace(const std::string& s) {

  const auto not_space = [](unsigned char c)->bool { return !std::isspace(c); };
  const auto first = std::find_if(s.begin(), s.end(), not_space);
  return {first, s.end()};

}

std::string kel::Utility::findAndReplaceAll(const std::string& source,
                                            const std::string& search,
                                            const std::string& replace) {

  if (search.empty()) {
    return source;
  }

  std::string modified;
  modified.reserve(source.size());
  std::size_t pos = 0;

  while (pos < source.size()) {

    const auto found = source.find(search, pos);
    if (found == std::string::npos) {
      modified.append(source, pos, std::string::npos);
      break;
    }
    modified.append(source, pos, found - pos);
    modified += replace;
    pos = found + search.size();
  }

  return modified;

}

// ---------------------------------------------------------------------------
// Tokenizers
// ---------------------------------------------------------------------------

std::vector<std::string_view> kel::Utility::viewTokenizer(const std::string_view& str_view, char delim) {

  std::vector<std::string_view> token_vector;
  std::size_t token_index = 0;

  for (std::size_t index = 0; index < str_view.size(); ++index) {

    if (str_view[index] == delim) {

      token_vector.emplace_back(str_view.substr(token_index, index - token_index));
      token_index = index + 1;

    }
  }

  token_vector.emplace_back(str_view.substr(token_index));
  return token_vector;

}

std::vector<std::string> kel::Utility::charTokenizer(const std::string_view& str_view, char delim) {

  const std::vector<std::string_view> views = viewTokenizer(str_view, delim);
  std::vector<std::string> str_vector;
  str_vector.reserve(views.size());

  for (std::string_view view : views) {

    str_vector.emplace_back(view);

  }
  return str_vector;

}

std::pair<std::string, std::string> kel::Utility::firstSplitChar(const std::string_view& source,
                                                                  bool(*char_delim_fn)(char c)) {

  const auto it = std::find_if(source.begin(), source.end(), char_delim_fn);
  return {std::string(source.begin(), it), std::string(it, source.end())};

}

// ---------------------------------------------------------------------------
// /proc parsing helpers
// ---------------------------------------------------------------------------

namespace {

// /proc/[pid]/stat format: pid (comm) state ppid ...
// 'comm' is parenthesized and may contain spaces, so naive operator>> parsing
// breaks. We read the whole line and locate the last ')' before parsing numbers.
bool read_proc_self_stat_line(std::ifstream& in, std::string& line) {
    return static_cast<bool>(std::getline(in, line)) && line.find_last_of(')') != std::string::npos;
}

} // anonymous namespace

// Returns virtual memory and resident set size in **GB**.
std::pair<double, double> kel::Utility::process_mem_usage() {

  static constexpr const char* STAT_STREAM = "/proc/self/stat";

  std::ifstream stat_stream(STAT_STREAM);
  if (!stat_stream.is_open()) {

    ExecEnv::log().error("Utility::process_mem_usage; Unable to open file stats: {}", STAT_STREAM);
    return {0.0, 0.0};

  }

  std::string line;
  if (!read_proc_self_stat_line(stat_stream, line)) {

    ExecEnv::log().error("Utility::process_mem_usage; Unable to read/parse file stats: {}", STAT_STREAM);
    return {0.0, 0.0};

  }

  std::istringstream rest(line.substr(line.find_last_of(')') + 1));

  std::string state;
  long ppid, pgrp, session, tty_nr, tpgid;
  unsigned long flags, minflt, cminflt, majflt, cmajflt;
  unsigned long utime, stime, cutime, cstime;
  long priority, nice;
  long num_threads, itrealvalue;
  unsigned long starttime;
  unsigned long vsize;
  long rss;

  rest >> state >> ppid >> pgrp >> session >> tty_nr
       >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
       >> utime >> stime >> cutime >> cstime >> priority >> nice
       >> num_threads >> itrealvalue >> starttime >> vsize >> rss;

  if (rest.fail()) {

    ExecEnv::log().error("Utility::process_mem_usage; Failed to parse stat fields");
    return {0.0, 0.0};

  }

  const long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024;
  double vm_usage = static_cast<double>(vsize) / 1024.0;
  double resident_set = static_cast<double>(rss) * page_size_kb;

    // Convert from kB to GB.
  vm_usage /= 1024.0 * 1024.0;
  resident_set /= 1024.0 * 1024.0;

  return {vm_usage, resident_set};

}

// Returns system and user CPU time in seconds.
std::pair<double, double> kel::Utility::process_time_usage() {

  static constexpr const char* STAT_STREAM = "/proc/self/stat";

  std::ifstream stat_stream(STAT_STREAM);
  if (!stat_stream.is_open()) {

    ExecEnv::log().error("Utility::process_time_usage; Unable to open file stats: {}", STAT_STREAM);
    return {0.0, 0.0};

  }

  std::string line;
  if (!read_proc_self_stat_line(stat_stream, line)) {

    ExecEnv::log().error("Utility::process_time_usage; Unable to read/parse file stats: {}", STAT_STREAM);
    return {0.0, 0.0};

  }

  std::istringstream rest(line.substr(line.find_last_of(')') + 1));

  std::string state;
  long ppid, pgrp, session, tty_nr, tpgid;
  unsigned long flags, minflt, cminflt, majflt, cmajflt;
  unsigned long utime, stime;

  rest >> state >> ppid >> pgrp >> session >> tty_nr
       >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
       >> utime >> stime;

  if (rest.fail()) {

    ExecEnv::log().error("Utility::process_time_usage; Failed to parse stat fields");
    return {0.0, 0.0};

  }

  return {static_cast<double>(stime) / 100.0, static_cast<double>(utime) / 100.0};

}

// Returns VM and physical memory in **MB** (matching /proc/status kB units).
std::pair<double, double> kel::Utility::process_mem_usage2() {

  static constexpr const char* STAT_STREAM = "/proc/self/status";
  static constexpr const char* VM_USAGE = "VmSize:";
  static constexpr const char* PM_USAGE = "VmRSS:";
  static constexpr double SCALE_MB = 1024.0;

  std::ifstream stat_stream(STAT_STREAM);

  if (!stat_stream.is_open()) {

    ExecEnv::log().error("Utility::process_mem_usage2; Unable to open file stats: {}", STAT_STREAM);
    return {0.0, 0.0};

  }

  const auto parse_kb_value = [](const std::string& line, const std::string& label) -> std::optional<double> {

    const auto label_pos = line.find(label);
    if (label_pos == std::string::npos) {
      return std::nullopt;
    }

    const auto digit_pos = line.find_first_of("0123456789", label_pos + label.size());
    if (digit_pos == std::string::npos) {
            return std::nullopt;
    }

    const auto end_pos = line.find_first_not_of("0123456789", digit_pos);
    const std::string number_str = line.substr(digit_pos, end_pos - digit_pos);

    try {

      return std::stod(number_str);

    } catch (const std::exception& e) {

      ExecEnv::log().error("process_mem_usage2; cannot parse {} with value: {} ({})", label, line, e.what());
      return std::nullopt;
    }

  };

  std::pair<double, double> mem_stats{0.0, 0.0};
  std::string stat_line;
  while (std::getline(stat_stream, stat_line)) {

    if (stat_line.find(VM_USAGE) != std::string::npos) {

      if (auto val = parse_kb_value(stat_line, VM_USAGE)) {

        mem_stats.first = *val;

      }

    } else if (stat_line.find(PM_USAGE) != std::string::npos) {

      if (auto val = parse_kb_value(stat_line, PM_USAGE)) {

        mem_stats.second = *val;

      }

    }

  }

  return {mem_stats.first / SCALE_MB, mem_stats.second / SCALE_MB};

}

// ---------------------------------------------------------------------------
// Statistics
// ---------------------------------------------------------------------------

std::pair<double, double> kel::Utility::stddev(const std::vector<double>& vec) {

  if (vec.empty()) {

    return {0.0, 0.0};

  }

  if (vec.size() == 1) {

    return {vec.front(), 0.0};

  }

  const double size = static_cast<double>(vec.size());
  const double mean = std::accumulate(vec.begin(), vec.end(), 0.0) / size;

  double sq_diff_sum = 0.0;

  for (const double val : vec) {

    const double diff = val - mean;
    sq_diff_sum += diff * diff;

  }

  const double variance = sq_diff_sum / (size - 1.0);
  return {mean, std::sqrt(variance)};

}

#else

// Given a file name as "path/file", returns "path"
std::string kel::Utility::filePath(const std::string& file_name) {

  fs::path file_path = fs::path(file_name);
  fs::path directory = file_path.parent_path();
  return directory.string();

}


// Returns the filename with the path directory appended to it "path/file", or appended subdirectory "path/sub_dir".
std::string kel::Utility::filePath(const std::string& file_name, const std::string& path) {

  fs::path directory_path = fs::path(path);
  fs::path file_path = directory_path / fs::path(file_name);
  return file_path.string();

}

// Utility to return a string with the current path plus file name, "path/file".
std::string kel::Utility::showPath() {

  return fs::current_path().string();

}

// Utility to append subdirectory "path/sub_dir" (same as above).
std::string kel::Utility::appendPath(const std::string& sub_directory, const std::string& path) { return filePath(sub_directory, path); }

// Check that a file exists at the file path.
bool kel::Utility::fileExists(const std::string& file_path) {

  bool file_exists = fs::is_regular_file(fs::path(file_path));
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

// If directory exists recreate it, else just create it.
bool kel::Utility::directoryRenew(const std::string& path) {

  if (directoryExists(path)) {

    return recreateDirectory(path);

  } else {

    return createDirectory(path);

  }

}


// Check that the directory exists at the specified path
bool kel::Utility::directoryExists(const std::string& path) {

  return fs::is_directory(fs::path(path));

}

// Create directory at the specified path, returns true if directory already exists.
bool kel::Utility::createDirectory(const std::string& path) {

  if (directoryExists(path)) {

    return true;

  }
  return fs::create_directory(fs::path(path));

}

// Recursively delete the contents of a directory and then the directory.
bool kel::Utility::deleteDirectory(const std::string& path_string) {

  fs::path fs_path(path_string);
  if (directoryExists(fs_path)) {

    return fs::remove_all(fs_path);

  } else {

    return false;

  }

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


// Only trim whitespace at the beginning of the string.
std::string kel::Utility::trimLeadingWhiteSpace(const std::string &s) {

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

  return start_trimmed_string;

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
std::pair<std::string, std::string> kel::Utility::firstSplitChar(const std::string_view& source, bool(* char_delim_fn)(char c)) {

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

// System and user CPU usage.
std::pair<double, double> kel::Utility::process_time_usage()
{

  // 'file' stat seems to give the most reliable results
  static const char* STAT_STREAM = "/proc/self/stat";


  std::ifstream stat_stream(STAT_STREAM,std::ios_base::in);

  if (not stat_stream.good()) {

    ExecEnv::log().error("Utility::process_time_usage; Unable to open file stats: {}", STAT_STREAM);
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

  double system_cpu = std::stod(stime) / 100.0;
  double user_cpu = std::stod(utime) / 100.0;

  return {system_cpu, user_cpu};

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

#endif
