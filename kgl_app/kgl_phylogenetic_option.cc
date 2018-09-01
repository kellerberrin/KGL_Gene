//
// Created by kellerberrin on 28/08/18.
//

#include "kgl_phylogenetic_option.h"
#include "kgl_exec_env.h"
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <sstream>


namespace kgl = kellerberrin::genome;
namespace bt = boost;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Reads and validates all the options defined in the specified option file.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::RuntimeOptions::getRuntimeOptionArray(const RuntimeOptionTag& option, RuntimeOptionVector& option_values) const {


  auto result = getMap().find(option);

  if (result == getMap().end()) {

    option_values.clear();
    return false;

  }

  option_values = result->second;
  return true;

}


bool kgl::RuntimeOptions::getRuntimeOption(const RuntimeOptionTag& option, std::string& option_value) const {

  RuntimeOptionVector option_values;
  bool result = getRuntimeOptionArray(option, option_values);

  if (not result) {

    return false;

  }

  if (option_values.empty()) {

    std::stringstream ss;
    printHelp(ss);
    ExecEnv::log().info("{}", ss.str());
    ExecEnv::log().critical("RuntimeOptions::getRuntimeOption; no argument available for option: {}", option);
    return false;

  }

  option_value = option_values.front();

  return true;

}




bool kgl::RuntimeOptions::readRuntimeOptions(const std::string& file_name, const std::string& work_directory) {

  // concatonate work directory and file number.

  std::string file_path = kgl::Utility::filePath(file_name, work_directory);

  std::ifstream readfile(file_path);

  if (not readfile.good()) {

    ExecEnv::log().critical("RuntimeOptions::readRuntimeOptions; could not open options file: {}", file_path);
    return false; // never reached.

  }

  bt::char_separator<char> item_key_sep(FIELD_DELIMITER_);
  size_t line_count = 0;
  size_t option_count = 0;

  while (true) {

    std::string record_str;
    if (std::getline(readfile, record_str).eof()) break;

    ++line_count;

    if ((record_str)[0] == COMMENT_CHAR_) continue;   // Ignore comment lines.

    std::vector<std::string> field_vector;
    bt::tokenizer<bt::char_separator<char>> tokenize_item(record_str, item_key_sep);

    for(auto iter_item = tokenize_item.begin(); iter_item != tokenize_item.end(); ++iter_item) {

      field_vector.push_back(*iter_item);

    }

    if (field_vector.empty()) continue;  // Skip any blank lines.

    ++option_count;

    RuntimeOptionTag option = field_vector.front();
    RuntimeOptionVector arguments = field_vector;
    arguments.erase(arguments.begin()); // Erase the first element (the option name).

    if (not parseOption(option, arguments, work_directory)) {

      ExecEnv::log().critical("RuntimeOptions::readRuntimeOptions; Error Parsing Options file: {}, line number: {} line: {}",
                              file_path, line_count, record_str);
      return false; // never reached.

    }

    std::pair<RuntimeOptionTag, RuntimeOptionVector> insert_option(option, arguments);
    auto result = runtime_option_map_.insert(insert_option);

    if (not result.second) {

      std::stringstream ss;
      printHelp(ss);
      ExecEnv::log().info("{}", ss.str());
      ExecEnv::log().critical("RuntimeOptions::readRuntimeOptions; Duplicate Option, file: {}, line number: {} line: {}",
                              file_path, line_count, record_str);

    }

  }

  if (not checkRequiredOptions()) {

    ExecEnv::log().critical("RuntimeOptions::readRuntimeOptions; Required option missing from options file: {}", file_path);
    return false; // never reached

  }

  readfile.close();
  return true;

}


bool kgl::RuntimeOptions::parseOption(const RuntimeOptionTag& option, RuntimeOptionVector& arguments, const std::string& work_directory) {

  PredefinedOptionType defined_option;
  if (not getDefinedOption(option, defined_option)) {

    ExecEnv::log().error("RuntimeOptions::parseOption; option: {} is unknown (not defined)", option);
    std::stringstream ss;
    printHelp(ss);
    ExecEnv::log().info("{}", ss.str());
    return false;
  }

  switch(defined_option.option_array_type) {

    case OptionArrayType::ARRAY: {
      if (arguments.empty()) {
        ExecEnv::log().error("RuntimeOptions::parseOption; option: {} argument type is an ARRAY, but no arguments supplied", option);
        std::stringstream ss;
        printHelp(ss);
        ExecEnv::log().info("{}", ss.str());
        return false;
      }
    }
      break;

    case OptionArrayType::SINGULAR: {
      if (arguments.size() != 1) {
        ExecEnv::log().error("RuntimeOptions::parseOption; option: {} argument type is SINGULAR. However, {} arguments supplied (1 expected).",
                             option, arguments.size());
        std::stringstream ss;
        printHelp(ss);
        ExecEnv::log().info("{}", ss.str());
        return false;
      }
    }
      break;

    case OptionArrayType::FLAG: {
      if (not arguments.empty()) {
        ExecEnv::log().error("RuntimeOptions::parseOption; option: {} argument type is FLAG. However, {} arguments supplied (zero expected).",
                             option, arguments.size());
        std::stringstream ss;
        printHelp(ss);
        ExecEnv::log().info("{}", ss.str());
        return false;
      }
    }
      break;

  }

  switch(defined_option.option_arument_type) {

    case OptionArgumentType::FILE_NAME: {

      for (auto& arg : arguments) {

        // If the arg is a filename then pre-pend the work directory to it.
        arg = kgl::Utility::filePath(arg, work_directory);

      }

    }

      break;

    case OptionArgumentType::STRING:
      break;

    case OptionArgumentType::INTEGER: {
      try {

        for (auto arg : arguments) {

          std::stoll(arg);

        }

      } catch (...) {

        ExecEnv::log().error("RuntimeOptions::parseOption; option: {} argument type is INTEGER. Bad argument type supplied.", option);
        std::stringstream ss;
        printHelp(ss);
        ExecEnv::log().info("{}", ss.str());
        return false;

      }
    }
      break;

    case OptionArgumentType::FLOAT: {
      try {

        for (auto arg : arguments) {

          std::stof(arg);

        }

      } catch (...) {

        ExecEnv::log().error("RuntimeOptions::parseOption; option: {} argument type is FLOAT. Bad argument type supplied.", option);
        std::stringstream ss;
        printHelp(ss);
        ExecEnv::log().info("{}", ss.str());
        return false;

      }
    }
      break;

  }

  return true;

}


