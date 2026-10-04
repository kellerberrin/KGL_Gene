//
// Created by kellerberrin on 15/2/20.
//

#include "kel_property_tree.h"
#include "kel_utility.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <charconv>
#include <sstream>
#include <fstream>
#include <iostream>
#include <shared_mutex>
#include <map>

namespace kel = kellerberrin;
namespace pt = boost::property_tree;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PropertyTree::PropertyImpl does all the heavy lifting using boost::property_tree.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using ImplSubTree = std::pair<std::string, kel::PropertyTree::PropertyImpl>;

/// The boost::property_tree PIMPL implementation. All accessor methods are const and do not
/// modify shared state, so a shared tree may be read concurrently from multiple threads.
class kel::PropertyTree::PropertyImpl {

public:


  PropertyImpl() = default;
  PropertyImpl(const pt::ptree& property_tree) : property_tree_(property_tree) {}
  PropertyImpl(const PropertyImpl&) =default;
  ~PropertyImpl() = default;

  [[nodiscard]] bool readPropertiesFile( const std::string& properties_file,
                                         const std::string& options_write_file,
                                         const std::string& parsed_write_file);

  [[nodiscard]] static bool writeProperties(const std::string& properties, const std::string& properties_write_file);

  [[nodiscard]] bool getProperty(const std::string& property_name, std::string& property) const;

  [[nodiscard]] bool getPropertyVector(const std::string& property_name, std::vector<std::string>& property_vector) const;

  [[nodiscard]] bool getNodeVector(const std::string& node_name, std::vector<std::string>& node_data_vector) const;

  [[nodiscard]] bool getProperty(const std::string& property_name, size_t& property) const;

  [[nodiscard]] std::stringstream treeTraversal() const;

  [[nodiscard]] bool checkProperty(const std::string& property_name) const;

  [[nodiscard]] bool getTreeVector(const std::string& property_name, std::vector<ImplSubTree>& tree_vector) const;

  [[nodiscard]] bool getTreeVector(std::vector<std::pair<std::string, PropertyImpl>>& tree_vector) const;

  template<class T> [[nodiscard]] T getData() const { return property_tree_.get_value<T>(); }

private:

  constexpr static std::string INCLUDE_TOKEN_{"#include"};  // The include file directive.
  constexpr static std::string LINE_IGNORE_{"//"};  // First 2 chars indicates a comment (pre-processor ignores line)
  constexpr static char LOGICAL_VAR_{'$'};  // Logical variable $<logical_variable>$ <replacement_text>
  constexpr static char INCLUDE_FILE_QUOTE_{'\"'};

  // boost property tree object
  pt::ptree property_tree_;
  std::unordered_map<std::string, std::string> logical_map_;

  void recursivePrintTree(std::stringstream& ss, const pt::ptree& property_tree, const std::string& parent, size_t depth) const;
  // These functions are recursive and throw file exceptions which are caught in readPropertiesFile().
  std::stringstream preProcessPropertiesFile(const std::string& properties_file);
  std::stringstream preProcessLogicalComments(const std::string& xml_file_name);  // Strips out comments and substitutes logical variables.
  void readRecursive(std::stringstream& ss, const std::string& properties_file, size_t& file_count);

};

// These functions throw file exceptions.
std::stringstream kel::PropertyTree::PropertyImpl::preProcessPropertiesFile(const std::string& properties_file) {

  std::stringstream ss;
  size_t file_count{0};

  readRecursive(ss, properties_file, file_count);

  ExecEnv::log().info("Runtime definition XML files parsed: {}", file_count);

  return ss;

}

// Strips out comments and substitutes logical variables.
std::stringstream kel::PropertyTree::PropertyImpl::preProcessLogicalComments(const std::string& xml_file_name) {

  static const size_t ignore_size = std::string(LINE_IGNORE_).size();
  std::ifstream xml_file(xml_file_name);

  if (not xml_file.good()) {

    ExecEnv::log().error("PropertyImpl::preProcessLogicalComments; Unable to open runtime XML file: {}", xml_file_name);
    throw std::runtime_error(xml_file_name);

  }

  std::stringstream scan_substitution;
  // First pass finds any defined logical substitutions.
  std::string line;
  while(std::getline(xml_file, line)) {

    // If the first two characters are "//" then the line is ignored.
    std::string ignore_text = line.substr(0, ignore_size);
    if (ignore_text == LINE_IGNORE_) {

      continue;

    }
    // If the first character is '$' then the line defines a logical definition
    if (line[0] == LOGICAL_VAR_) {

      auto line_split = Utility::viewTokenizer(line, LOGICAL_VAR_);
      if (line_split.size() != 3) {

        ExecEnv::log().warn("$<logical_variable>$ \"<replacement_text>\"; Invalid logical format: {}, tokens: {}", line, line_split.size());
        continue;
      }

      auto logical_var = line_split[1];
      auto replacement_vec = Utility::viewTokenizer(line_split[2], INCLUDE_FILE_QUOTE_);
      if (replacement_vec.size() < 2) {

        ExecEnv::log().warn("\"<replacement_text>\"; Invalid format: {}", line_split[1]);
        continue;

      }
      auto replacement = replacement_vec[1];

      ExecEnv::log().info("Substitute: ${}$ -> \"{}\"", logical_var, replacement);
      logical_map_.emplace(std::string(logical_var),std::string(replacement));

      continue;

    }

    scan_substitution << line << '\n';

  }

  // Second pass substitutes any logical definitions with replacement text.
  std::stringstream substituted_xml;
  while(std::getline(scan_substitution, line)) {

    auto line_split = Utility::viewTokenizer(line, LOGICAL_VAR_);
    // A logical var will split the line into 3 tokens <text1>$<logical_var>$<text2>.
    if (line_split.size() == 3) {

      auto logical_iter = logical_map_.find(std::string(line_split[1]));
      // Check if not found
      if (logical_iter == logical_map_.end()) {

        ExecEnv::log().warn("Logical variable ${}$ not found for line: {} - no substitution performed", line_split[1], line);

      } else {

        // Reconstruct the substituted line.
        auto [logical_var, replacement] = *logical_iter;
        line = std::string(line_split[0]) + replacement + std::string(line_split[2]);

      }

    }

    substituted_xml << line << '\n';

  }

  return substituted_xml;

}

// This function preprocesses the runtime XML file by allowing a c++ style include directive, e.g. #include "subdir/include.xml".
// This allows the runtime XML file to be broken up and simplified.
// Redundant include statements can be disabled by prefixing with '//' in the first two characters of the line.
// For example '//#include "subdir/include.xml' is a disabled include statement.
void kel::PropertyTree::PropertyImpl::readRecursive(std::stringstream& ss, const std::string& properties_file, size_t& file_count) {

  static const size_t include_token_size = std::string(INCLUDE_TOKEN_).size();

  ++file_count;
  auto logicalProcessed = preProcessLogicalComments(properties_file);

  std::string line;
  while(std::getline(logicalProcessed, line))
  {

    if (line.contains(INCLUDE_TOKEN_)) {

      std::string trimmed_line = Utility::trimLeadingWhiteSpace(line);
      if (trimmed_line.starts_with(INCLUDE_TOKEN_)) {

        // The include XML file spec should be in quotes, e.g #include "subdir/include.xml".
        std::string file_spec = Utility::trimAllWhiteSpace(line.substr(include_token_size));
        file_spec = Utility::trimAllChar(file_spec, INCLUDE_FILE_QUOTE_);
        file_spec = Utility::filePath(file_spec, Utility::filePath(properties_file));
        // Recursively include XML file.
        readRecursive(ss, file_spec, file_count);


      } else {

        ExecEnv::log().warn("XML property file: {} contains malformed {} directive: {}", properties_file, INCLUDE_TOKEN_, line);
        ss << line << '\n';

      }

    } else {

      ss << line << '\n';

    }

  }

}


bool kel::PropertyTree::PropertyImpl::readPropertiesFile( const std::string& properties_file,
                                                          const std::string& options_write_file,
                                                          const std::string& parsed_write_file) {

  try {

    std::stringstream ss = preProcessPropertiesFile(properties_file);

    // Write the text xml file (prior to parsing) if the option_file_out argument is specified.
    if (not options_write_file.empty()) {

      if (not writeProperties(ss.str(), options_write_file)) {

        ExecEnv::log().error("Cannot write processed options file {}", options_write_file);

      }

    }

    // Parse the preprocessed text file into an xml tree.
    pt::read_xml(ss, property_tree_);

    // If an xml output file specified, then write out the parsed xml file.
    if (not parsed_write_file.empty()) {

      if (not writeProperties(treeTraversal().str(), parsed_write_file)) {

        ExecEnv::log().error("Cannot write parsed xml file {}", parsed_write_file);

      }

    }

  } catch(const std::exception& e) {

    ExecEnv::log().error("PropertyTree; Missing or Malformed property tree in file: {}, error: {}", properties_file, e.what());
    return false;

  }

  return true;

}

bool kel::PropertyTree::PropertyImpl::writeProperties(const std::string& properties, const std::string& properties_write_file) {

  std::ofstream properties_file(properties_write_file);

  if (not properties_file.good()) {

    ExecEnv::log().error("PropertyImpl::writePropertiesFile; could not open file: {} for properties file output", properties_write_file);
    return false;

  }

  properties_file << properties;

  return properties_file.good();

}


bool kel::PropertyTree::PropertyImpl::checkProperty(const std::string& property_name) const {

  try {

    property_tree_.get<std::string>(property_name);

  }
  catch (...) {

    return false;

  }

  return true;

}


bool kel::PropertyTree::PropertyImpl::getProperty(const std::string& property_name,  std::string& property) const {

  try {

    property = property_tree_.get<std::string>(property_name);
    property = Utility::trimEndWhiteSpace(property);
    if (property.empty()) {

      ExecEnv::log().error("PropertyTree; Well-formed Property Tree but Property: {} not found or NULL", property_name);
      return false;

    }

  }
  catch (const std::exception& e) {

    ExecEnv::log().error("Exception: PropertyTree::PropertyImpl::getProperty; Property: {} not found, error: {}", property_name, e.what());
    ExecEnv::log().error("*********** Property Tree Contents *************");
    ExecEnv::log().error(treeTraversal().str());
    ExecEnv::log().error("**********************************************");
    return false;

  }

  return true;

}


bool kel::PropertyTree::PropertyImpl::getPropertyVector(const std::string& property_name, std::vector<std::string>& property_vector) const {


  try {

    for (auto const& child : property_tree_.get_child(property_name)) {

      // The data function is used to access the data stored in a node.
      property_vector.push_back(Utility::trimEndWhiteSpace(child.second.data()));

    }

  }
  catch (const std::exception& e) {

    ExecEnv::log().error("PropertyTree::getPropertyVector; Property Vector: {} not found, error: {}", property_name, e.what());
    ExecEnv::log().error("***********Property Tree Contents*************");
    ExecEnv::log().error(treeTraversal().str());
    ExecEnv::log().error("**********************************************");
    return false;

  }

  return true;

}



bool kel::PropertyTree::PropertyImpl::getNodeVector(const std::string& node_name, std::vector<std::string>& node_vector) const {

  try {

    for (auto const& tree : property_tree_) {

      if (tree.first == node_name) {

        node_vector.push_back(tree.second.get_value<std::string>());

      }

    }

  }
  catch (const std::exception& e) {

    ExecEnv::log().error("PropertyTree::getPropertyVector; Property Vector: {} not found, error: {}", node_name, e.what());
    ExecEnv::log().error("***********Property Tree Contents*************");
    ExecEnv::log().error(treeTraversal().str());
    ExecEnv::log().error("**********************************************");
    return false;

  }

  return true;

}

bool kel::PropertyTree::PropertyImpl::getTreeVector(std::vector<std::pair<std::string, PropertyImpl>>& tree_vector) const {

  tree_vector.clear();

  try {

    for (auto const& sub_tree : property_tree_) {

      tree_vector.emplace_back(ImplSubTree(sub_tree.first, PropertyImpl(sub_tree.second)));

    }

  }
  catch (...) {

    // No sub-tree is not an error.
    return true;

  }

  return true;

}




bool kel::PropertyTree::PropertyImpl::getTreeVector(const std::string& property_name, std::vector<std::pair<std::string, PropertyImpl>>& tree_vector) const {

  tree_vector.clear();

  try {

    for (auto const& sub_tree : property_tree_.get_child(property_name)) {

      tree_vector.emplace_back(ImplSubTree(sub_tree.first, PropertyImpl(sub_tree.second)));

    }

  }
  catch (...) {

    // No sub-tree is not an error.
    return true;

  }

  return true;

}


bool kel::PropertyTree::PropertyImpl::getProperty(const std::string& property_name, size_t& property) const {

  std::string size_string;
  if (not getProperty(property_name, size_string)) {

    return false;

  }

  // Use std::from_chars for no-throw numeric parsing and full-consumption validation.
  size_t parsed_value{0};
  const auto [ptr, ec] = std::from_chars(size_string.data(), size_string.data() + size_string.size(), parsed_value);
  if (ec != std::errc() or ptr != size_string.data() + size_string.size()) {

    ExecEnv::log().error("PropertyTree; Property: {}, Value: {} is not an unsigned integer", property_name, size_string);
    return false;

  }

  property = parsed_value;
  return true;

}


std::stringstream kel::PropertyTree::PropertyImpl::treeTraversal() const {

  std::stringstream ss;
  recursivePrintTree(ss, property_tree_, "", 0);

  return ss;

}



void kel::PropertyTree::PropertyImpl::recursivePrintTree( std::stringstream& ss,
                                                          const pt::ptree& property_tree,
                                                          const std::string& parent,
                                                          size_t depth) const {

  ++depth; // Next tree branch.
  const size_t horiz_spaces = 2;  // Spaces per tree branch.
  size_t spaces = horiz_spaces * depth;

  for (const auto& item : property_tree) {

    std::string key = item.first;
    std::string parent_key;
    if (parent.empty()) {

      parent_key = key;

    } else {

      parent_key = parent + "." + key;

    }
    std::string value = item.second.data();
    value = Utility::trimEndWhiteSpace(value);
    if (not value.empty()) {

      ss << std::format("{: <{}}{}={}\n", "", spaces, parent_key, value);

    } else {

      ss << std::format("{: <{}}{}\n", "", spaces, parent_key);

    }
    // Recursive call to descend through the xml tree.
    recursivePrintTree(ss, item.second, parent_key, depth);

    ss << std::format("{: <{}}/{}\n", "", spaces, parent_key);

  }

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PropertyTree is a public facade class that passes the functionality onto PropertyTree::PropertyImpl.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kel::PropertyTree::PropertyTree() : properties_impl_ptr_(std::make_unique<kel::PropertyTree::PropertyImpl>()) {}

kel::PropertyTree::PropertyTree(const PropertyTree& property_tree) {

  std::shared_lock lock(property_tree.tree_mutex_);
  properties_impl_ptr_ = std::make_unique<kel::PropertyTree::PropertyImpl>(*(property_tree.properties_impl_ptr_));

}

kel::PropertyTree::PropertyTree(const PropertyImpl& properties_impl) : properties_impl_ptr_(std::make_unique<kel::PropertyTree::PropertyImpl>(properties_impl)) {}

kel::PropertyTree::~PropertyTree() {}  // Required because of incomplete pimpl type.

// Functionality passed to the implmentation.

bool kel::PropertyTree::readProperties( const std::string& properties_file,
                                        const std::string& options_write_file,
                                        const std::string& parsed_write_file) {

  auto new_impl = std::make_unique<PropertyImpl>();
  if (not new_impl->readPropertiesFile(properties_file, options_write_file, parsed_write_file)) {

    return false;

  }

  std::unique_lock lock(tree_mutex_);
  properties_impl_ptr_ = std::move(new_impl);
  return true;

}


bool kel::PropertyTree::getProperty(const std::string& property_name, std::string& property) const {

  std::shared_lock lock(tree_mutex_);
  return properties_impl_ptr_->getProperty(property_name, property);

}


bool kel::PropertyTree::getProperty(const std::string& property_name, size_t& property) const {

  std::shared_lock lock(tree_mutex_);
  return properties_impl_ptr_->getProperty(property_name, property);

}


bool kel::PropertyTree::getOptionalProperty(const std::string& property_name, std::string& property) const {

  std::shared_lock lock(tree_mutex_);
  if (properties_impl_ptr_->checkProperty(property_name)) {

    return properties_impl_ptr_->getProperty(property_name, property);

  } else {

    return false;

  }

}

std::stringstream kel::PropertyTree::treeTraversal() const {

  std::shared_lock lock(tree_mutex_);
  return properties_impl_ptr_->treeTraversal();

}


bool kel::PropertyTree::checkProperty(const std::string& property_name) const {

  std::shared_lock lock(tree_mutex_);
  return properties_impl_ptr_->checkProperty(property_name);

}


bool kel::PropertyTree::getFileProperty(const std::string& property_name, const std::string& work_directory, std::string& file_path) const {

  std::shared_lock lock(tree_mutex_);
  std::string file_name;

  if (not properties_impl_ptr_->getProperty(property_name, file_name)) {

    ExecEnv::log().warn("PropertyTree::getFileProperty; Requested file property: {} not found. A list of all valid properties follows:",
                        property_name);

    auto ss = properties_impl_ptr_->treeTraversal();
    std::cout << ss.str() << std::endl;

    return false;

  }

  file_path = Utility::filePath(file_name, work_directory);

  if (not Utility::fileExists(file_path)) {

    ExecEnv::log().warn("PropertyTree::getFileProperty; File: {} does not exist", file_path);
    return false;

  }

  return true;

}



bool kel::PropertyTree::getFileCreateProperty(const std::string& property_name, const std::string& work_directory, std::string& file_path) const {

  std::shared_lock lock(tree_mutex_);
  std::string file_name;

  if (not properties_impl_ptr_->getProperty(property_name, file_name)) {

    ExecEnv::log().warn("PropertyTree::getFileCreateProperty; Requested file property: {} not found. A list of all valid properties follows:",
                        property_name);

    auto ss = properties_impl_ptr_->treeTraversal();
    std::cout <<ss.str() << std::endl;

    return false;

  }

  file_path = Utility::filePath(file_name, work_directory);

  if (Utility::fileExists(file_path)) {

    return true;

  }

  if (Utility::fileExistsCreate(file_path)) {

    ExecEnv::log().info("Requested Property; Created file: {}", file_path);
    return true;

  } else {

    ExecEnv::log().warn("PropertyTree::getFileCreateProperty; File: {} does not exist and could not be created", file_path);
    return false;

  }

}



bool kel::PropertyTree::getOptionalFileProperty(const std::string& property_name, const std::string& work_directory, std::string& file_path) const {

  std::shared_lock lock(tree_mutex_);
  std::string file_name;

  if (not (properties_impl_ptr_->checkProperty(property_name)
           and properties_impl_ptr_->getProperty(property_name, file_name))) {

    return false;

  }

  file_path = Utility::filePath(file_name, work_directory);

  if (not Utility::fileExists(file_path)) {

    ExecEnv::log().critical("PropertyTree::getFileProperty; Specified optional File: {} does not exist, remove from XML specification.", file_path);
    return false;

  }

  return true;

}


bool kel::PropertyTree::getPropertyVector(const std::string& property_name, std::vector<std::string>& property_vector) const {

  std::shared_lock lock(tree_mutex_);
  return properties_impl_ptr_->getPropertyVector(property_name, property_vector);

}


bool kel::PropertyTree::getNodeVector(const std::string& node_name, std::vector<std::string>& node_vector) const {

  std::shared_lock lock(tree_mutex_);
  return properties_impl_ptr_->getNodeVector(node_name, node_vector);

}


bool kel::PropertyTree::getPropertyTreeVector(const std::string& property_name, std::vector<SubPropertyTree>& property_tree_vector) const {

  property_tree_vector.clear();

  std::vector<ImplSubTree> tree_vector;
  {
    std::shared_lock lock(tree_mutex_);
    if (not properties_impl_ptr_->getTreeVector(property_name, tree_vector)) {

      return false;

    }
  }

  for (auto const& [sub_tree_tag, sub_tree] : tree_vector) {

    // Ignore comment  and help tags.
    if (sub_tree_tag != COMMENT_ and sub_tree_tag != HELP_) {

      property_tree_vector.emplace_back(SubPropertyTree(sub_tree_tag, PropertyTree(sub_tree)));

    }

  }

  return true;

}


bool kel::PropertyTree::getPropertySubTreeVector(std::vector<SubPropertyTree>& property_tree_vector) const {

  property_tree_vector.clear();

  std::vector<ImplSubTree> tree_vector;
  {
    std::shared_lock lock(tree_mutex_);
    if (not properties_impl_ptr_->getTreeVector(tree_vector)) {

      return false;

    }
  }

  for (auto const& [property_tree_id, property_tree] : tree_vector) {

    property_tree_vector.emplace_back(SubPropertyTree(property_tree_id, PropertyTree(property_tree)));

  }

  return true;

}



std::string kel::PropertyTree::getValue() const {

  std::shared_lock lock(tree_mutex_);
  try {

    return properties_impl_ptr_->getData<std::string>();

  }
  catch (const std::exception& e) {

    ExecEnv::log().error("PropertyTree::getValue; error: {}", e.what());
    return "";

  }
  catch (...) {

    ExecEnv::log().error("PropertyTree::getValue; unknown error");
    return "";

  }

}
