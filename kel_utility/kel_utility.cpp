// Copyright 2023 Kellerberrin
//


#include "kel_utility.h"
#include "kel_exec_env.h"

#include <algorithm>
#include <charconv>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
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

// ---------------------------------------------------------------------------
// Path / filesystem helpers
// ---------------------------------------------------------------------------

/// Input is "path/file.ext" function returns "path".
std::string kel::Utility::filePath(const std::string& file_name) {

  return fs::path(file_name).parent_path().string();

}

/// Utility to concatenate "path/file.ext".
std::string kel::Utility::filePath(const std::string& file_name, const std::string& path) {

  return (fs::path(path) / file_name).string();

}

/// Utility to return a string with the current path plus file name.
std::string kel::Utility::showPath() {

  std::error_code ec;
  auto path = fs::current_path(ec);

  if (ec) {

    ExecEnv::log().error("Utility::showPath; error getting current path: {}", ec.message());
    return {};

  }

  return path.string();

}

/// Utility to append subdirectory "path/sub_dir" (same as above).
std::string kel::Utility::appendPath(const std::string& sub_directory, const std::string& path) {

  return filePath(sub_directory, path);

}

/// Check that a file exists at the file path.
bool kel::Utility::fileExists(const std::string& file_path) {

  std::error_code ec;
  bool exists = fs::is_regular_file(file_path, ec);

  if (ec) {

    ExecEnv::log().error("Utility::fileExists; error checking {}: {}", file_path, ec.message());
    return false;

  }
  return exists;

}

/// Check that a file exists at the file path, creates zero sized file if not.
bool kel::Utility::fileExistsCreate(const std::string& file_path) {

  if (fileExists(file_path)) {

    return true;

  }
  std::ofstream empty_file(file_path);
  if (not empty_file.is_open()) {

    ExecEnv::log().error("Utility::fileExistsCreate; failed to create {}", file_path);
    return false;

  }
  return true;

}

/// Check that the directory exists at the specified path.
bool kel::Utility::directoryExists(const std::string& path) {

  std::error_code ec;
  bool exists = fs::is_directory(fs::path(path), ec);

  if (ec) {

    ExecEnv::log().error("Utility::directoryExists; error checking {}: {}", path, ec.message());
    return false;

  }
  return exists;

}

/// Create directory at the specified path.
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

/// Recursively deletes the contents of a directory and all subdirectories.
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

/// Delete the directory and its contents and then recreate the directory.
bool kel::Utility::recreateDirectory(const std::string& path) {

  if (!deleteDirectory(path)) {

    return false;
  }

  return createDirectory(path);

}

/// If directory exists, recreate it, else just create it.
bool kel::Utility::directoryRenew(const std::string& path) {

  if (directoryExists(path)) {

    return recreateDirectory(path);

  }

  return createDirectory(path);

}

/// Input is "path/file.ext", returns "ext".
std::string kel::Utility::fileExtension(const std::string& file_name) {

  return fs::path(file_name).extension().string();

}

/// Input is path/file.ext", returns "file".
std::string kel::Utility::fileName(const std::string& file_name) {

  return fs::path(file_name).filename().string();

}

// ---------------------------------------------------------------------------
// Environment / string helpers
// ---------------------------------------------------------------------------

/// Translate a linux environment variable.
std::optional<std::string> kel::Utility::getEnvironment(const std::string& env_var) {
    if (const char* val = std::getenv(env_var.c_str())) {
        return std::string(val);
    }
    return std::nullopt;
}

/// Covert to upper case.
std::string kel::Utility::toupper(const std::string& s) {

  std::string upper_string;
  upper_string.reserve(s.size());

  auto toupper_lambda = [](unsigned char c)->char { return static_cast<char>(std::toupper(c)); };
  std::transform(s.begin(), s.end(), std::back_inserter(upper_string), toupper_lambda);

  return upper_string;
}

/// Trim any whitespace in a string.
std::string kel::Utility::trimAllWhiteSpace(const std::string& s) {

  std::string clean_string;
  clean_string.reserve(s.size());

  auto ws_pred = [](unsigned char c)->bool { return !std::isspace(c); };
  std::copy_if(s.begin(), s.end(), std::back_inserter(clean_string), ws_pred);
  return clean_string;

}

/// Returns a string with all nc char removed.
std::string kel::Utility::trimAllChar(const std::string& s, char nc) {

  std::string mod_string;
  mod_string.reserve(s.size());

  const auto target = static_cast<unsigned char>(nc);
  auto trim_pred = [target](unsigned char c)->bool { return c != target; };
  std::copy_if(s.begin(), s.end(), std::back_inserter(mod_string), trim_pred);

  return mod_string;

}

/// Only trim whitespace at either end of the string.
std::string kel::Utility::trimEndWhiteSpace(const std::string& s) {

  const auto not_space = [](unsigned char c)->bool { return !std::isspace(c); };
  const auto first = std::find_if(s.begin(), s.end(), not_space);
  if (first == s.end()) {
        return {};
  }

  const auto last = std::find_if(s.rbegin(), s.rend(), not_space).base();
  return {first, last};

}

/// Only trim whitespace at beginning of the string.
std::string kel::Utility::trimLeadingWhiteSpace(const std::string& s) {

  const auto not_space = [](unsigned char c)->bool { return !std::isspace(c); };
  const auto first = std::find_if(s.begin(), s.end(), not_space);
  return {first, s.end()};

}

/// Replace all occurrences of search with replace in source.
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

/// Tokenize a string using delimiter chars, return std::string_view tokens.
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

/// Tokenize a string using delimiter chars, return std::string tokens.
std::vector<std::string> kel::Utility::charTokenizer(const std::string_view& str_view, char delim) {

  const std::vector<std::string_view> views = viewTokenizer(str_view, delim);
  std::vector<std::string> str_vector;
  str_vector.reserve(views.size());

  for (std::string_view view : views) {

    str_vector.emplace_back(view);

  }
  return str_vector;

}

/// Split string on encountering a specific char. Default version splits on first whitespace.
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

/// Returns virtual memory and resident set size in **GB**.
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

/// Returns system and user CPU time in seconds.
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

/// Returns VM and physical memory in **MB** (matching /proc/status kB units).
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

    double value{0.0};
    const auto [ptr, ec] = std::from_chars(number_str.data(), number_str.data() + number_str.size(), value);
    if (ec != std::errc() || ptr != number_str.data() + number_str.size()) {

      ExecEnv::log().error("process_mem_usage2; cannot parse {} with value: {}", label, line);
      return std::nullopt;

    }

    return value;

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

/// mean = pair.first,  sample stdev = pair.second.
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
