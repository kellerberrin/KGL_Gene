//
// Created by kellerberrin on 15/2/20.
//

#include "kel_property_tree.h"
#include "kel_utility.h"

// To stop boost complaining about global placeholders on boost/bind.
#define BOOST_BIND_GLOBAL_PLACEHOLDERS 1
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/spirit/include/lex_lexertl.hpp>

#include <sstream>
#include <stack>

namespace kel = kellerberrin;
namespace pt = boost::property_tree;
namespace lex = boost::spirit::lex;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PropertyTree::PropertyImpl does all the heavy lifting using boost::property_tree.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using ImplSubTree = std::pair<std::string, kel::PropertyTree::PropertyImpl>;

class kel::PropertyTree::PropertyImpl {

public:


  PropertyImpl() = default;
  PropertyImpl(const pt::ptree& property_tree) : property_tree_(property_tree) {}
  PropertyImpl(const PropertyImpl&) =default;
  ~PropertyImpl() = default;

  bool readPropertiesFile(const std::string& properties_file);

//  bool getProperty(const std::string& property_name, std::string& property) const { return getProperty(property_tree_, property_name, property); }

  bool getProperty(const std::string& property_name, std::string& property) const;

  bool getPropertyVector(const std::string& property_name, std::vector<std::string>& property_vector) const;

  bool getNodeVector(const std::string& node_name, std::vector<std::string>& node_data_vector) const;

  bool getProperty(const std::string& property_name, size_t& property) const;


  void treeTraversal() const;

  bool checkProperty(const std::string& property_name) const;

  bool getTreeVector(const std::string& property_name, std::vector<ImplSubTree>& tree_vector) const;

  bool getTreeVector(std::vector<std::pair<std::string, PropertyImpl>>& tree_vector) const;

    template<class T> T getData() const { return property_tree_.get_value<T>(); }

private:

  constexpr static const char* JSON_FILE_EXT_{"JSON"};
  constexpr static const char* INCLUDE_TOKEN_{"#include"};
  constexpr static const char* LINE_IGNORE_{"//"};
  constexpr static const char INCLUDE_FILE_QUOTE_{'\"'};

  // boost property tree object
  pt::ptree property_tree_;

  void printTree(const std::string& parent, const pt::ptree& property_tree) const;
  bool getPropertyTree(const pt::ptree& property_tree, const std::string& property_name, std::string& property) const;
  // These functions are recursive and throw file exceptions which are caught in getPropertyTree().
  std::stringstream preProcessPropertiesFile(const std::string& properties_file);
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

// This function preprocesses the runtime XML file by allowing a c++ style include directive, e.g. #include "subdir/include.xml".
// This allows the runtime XML file to be broken up and simplified.
// Redundant include statements can be disabled by prefixing with '//' in the first two characters of the line.
// For example '//#include "subdir/include.xml' is a disabled include statement.
void kel::PropertyTree::PropertyImpl::readRecursive(std::stringstream& ss, const std::string& properties_file, size_t& file_count) {

  static const size_t include_token_size = std::string(INCLUDE_TOKEN_).size();
  static const size_t ignore_size = std::string(LINE_IGNORE_).size();

  ++file_count;
  std::ifstream xml_file(properties_file);

  if (not xml_file.good()) {

    ExecEnv::log().error("PropertyImpl::preProcessPropertiesFile; Unable to open runtime XML file: {}", properties_file);
    throw std::runtime_error(properties_file);

  }

  std::string line;
  while(std::getline(xml_file, line))
  {
    // If the first two characters are "//" then the line is ignored.
    std::string ignore_text = line.substr(0, ignore_size);
    if (ignore_text == LINE_IGNORE_) {

      continue;

    }

    // Recursively include XML files.
    if (line.find(INCLUDE_TOKEN_) != std::string::npos) {

      // The include XML file spec should be in quotes, e.g #include "subdir/include.xml".
      std::string file_spec = Utility::trimAllWhiteSpace(line.substr(include_token_size));
      file_spec = Utility::trimAllChar(file_spec, INCLUDE_FILE_QUOTE_);
      file_spec = Utility::filePath(file_spec, Utility::filePath(properties_file));
      readRecursive(ss, file_spec, file_count);

    } else {

      ss << line << '\n';

    }

  }

}


bool kel::PropertyTree::PropertyImpl::readPropertiesFile(const std::string& properties_file) {

  std::string uc_file_extension = Utility::toupper(Utility::fileExtension(properties_file));

  try {

    if (uc_file_extension == JSON_FILE_EXT_) {

      pt::read_json(properties_file, property_tree_);

    } else {

      std::stringstream ss = preProcessPropertiesFile(properties_file);
      pt::read_xml(ss, property_tree_);

    }

  }
  catch(std::runtime_error& e) {

    ExecEnv::log().error("PropertyTree; Missing or Malformed property tree in file: {}, error: {}", properties_file, e.what());
    return false;

  }

  return true;

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
  catch (const std::runtime_error& e) {

    ExecEnv::log().error("Exception: PropertyTree::PropertyImpl::getProperty; Property: {} not found, error: {}", property_name, e.what());
    ExecEnv::log().error("*********** Property Tree Contents *************");
    treeTraversal();
    ExecEnv::log().error("**********************************************");
    return false;

  }

  return true;

}


bool kel::PropertyTree::PropertyImpl::getPropertyTree(const pt::ptree& property_tree, const std::string& property_name,  std::string& property) const {

  try {

    property = property_tree.get<std::string>(property_name);
    property = Utility::trimEndWhiteSpace(property);

  }
  catch (const std::runtime_error& e) {

    ExecEnv::log().error("PropertyTree; Property: {} not found, error: {}", property_name, e.what());
    ExecEnv::log().error("***********Property Tree Contents*************");
    treeTraversal();
    ExecEnv::log().error("**********************************************");
    return false;

  }

  return true;

}

bool kel::PropertyTree::PropertyImpl::getPropertyVector(const std::string& property_name, std::vector<std::string>& property_vector) const {


  try {

    for (auto child : property_tree_.get_child(property_name)) {

      // The data function is used to access the data stored in a node.
      property_vector.push_back(Utility::trimEndWhiteSpace(child.second.data()));

    }

  }
  catch (const std::runtime_error& e) {

    ExecEnv::log().error("PropertyTree::getPropertyVector; Property Vector: {} not found, error: {}", property_name, e.what());
    ExecEnv::log().error("***********Property Tree Contents*************");
    treeTraversal();
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
  catch (const std::runtime_error& e) {

    ExecEnv::log().error("PropertyTree::getPropertyVector; Property Vector: {} not found, error: {}", node_name, e.what());
    ExecEnv::log().error("***********Property Tree Contents*************");
    treeTraversal();
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

    for (auto sub_tree : property_tree_.get_child(property_name)) {

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

  try {

    property= std::stoull(size_string);

  }
  catch (const std::runtime_error& e) {

    ExecEnv::log().error("PropertyTree; Property: {}, Value: {} is not an unsigned integer, error: {}", property_name, size_string, e.what());
    return false;

  }

  return true;

}


void kel::PropertyTree::PropertyImpl::treeTraversal() const {

  printTree("", property_tree_);

}



void kel::PropertyTree::PropertyImpl::printTree(const std::string& parent, const pt::ptree& property_tree) const {


  for (auto item : property_tree) {

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

      ExecEnv::log().info("PropertyTree; key: {}, value: {}", parent_key, value);

    } else {

      ExecEnv::log().info("PropertyTree; key: {}", parent_key);

    }
    printTree(parent_key, item.second);

  }

}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PropertyTree is a public facade class that passes the functionality onto PropertyTree::PropertyImpl.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kel::PropertyTree::PropertyTree() : properties_impl_ptr_(std::make_unique<kel::PropertyTree::PropertyImpl>()) {}

kel::PropertyTree::PropertyTree(const PropertyTree& property_tree) : properties_impl_ptr_(std::make_unique<kel::PropertyTree::PropertyImpl>(*(property_tree.properties_impl_ptr_))) {}

kel::PropertyTree::PropertyTree(const PropertyImpl& properties_impl) : properties_impl_ptr_(std::make_unique<kel::PropertyTree::PropertyImpl>(properties_impl)) {}

kel::PropertyTree::~PropertyTree() {}  // DO NOT DELETE or USE DEFAULT. Required because of incomplete pimpl type.

// Functionality passed to the implmentation.

bool kel::PropertyTree::readProperties(const std::string& properties_file) {

  return properties_impl_ptr_->readPropertiesFile(properties_file);

}


bool kel::PropertyTree::getProperty(const std::string& property_name, std::string& property) const {

  return properties_impl_ptr_->getProperty(property_name, property);

}


bool kel::PropertyTree::getProperty(const std::string& property_name, size_t& property) const {

  return properties_impl_ptr_->getProperty(property_name, property);

}


bool kel::PropertyTree::getOptionalProperty(const std::string& property_name, std::string& property) const {

  if (checkProperty(property_name)) {

    return properties_impl_ptr_->getProperty(property_name, property);

  } else {

    return false;

  }

}

void kel::PropertyTree::treeTraversal() const {

  return properties_impl_ptr_->treeTraversal();

}


bool kel::PropertyTree::checkProperty(const std::string& property_name) const {

  return properties_impl_ptr_->checkProperty(property_name);

}


bool kel::PropertyTree::getFileProperty(const std::string& property_name, const std::string& work_directory, std::string& file_path) const {

  std::string file_name;

  if (not getProperty(property_name, file_name)) {

    ExecEnv::log().warn("PropertyTree::getFileProperty; Requested file property: {} not found. A list of all valid properties follows:",
                        property_name);

    treeTraversal();
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

  std::string file_name;

  if (not getProperty(property_name, file_name)) {

    ExecEnv::log().warn("PropertyTree::getFileCreateProperty; Requested file property: {} not found. A list of all valid properties follows:",
                        property_name);

    treeTraversal();
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

  std::string file_name;

  if (not getOptionalProperty(property_name, file_name)) {

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

  return properties_impl_ptr_->getPropertyVector(property_name, property_vector);

}


bool kel::PropertyTree::getNodeVector(const std::string& node_name, std::vector<std::string>& node_vector) const {

  return properties_impl_ptr_->getNodeVector(node_name, node_vector);

}


bool kel::PropertyTree::getPropertyTreeVector(const std::string& property_name, std::vector<SubPropertyTree>& property_tree_vector) const {

  property_tree_vector.clear();

  std::vector<ImplSubTree> tree_vector;
  if (not properties_impl_ptr_->getTreeVector(property_name, tree_vector)) {

    return false;

  }

  for (auto impl_tree : tree_vector) {

    property_tree_vector.emplace_back(SubPropertyTree(impl_tree.first, PropertyTree(impl_tree.second)));

  }

  return true;

}


bool kel::PropertyTree::getPropertySubTreeVector(std::vector<SubPropertyTree>& property_tree_vector) const {

  property_tree_vector.clear();

  std::vector<ImplSubTree> tree_vector;
  if (not properties_impl_ptr_->getTreeVector(tree_vector)) {

    return false;

  }

  for (auto impl_tree : tree_vector) {

    property_tree_vector.emplace_back(SubPropertyTree(impl_tree.first, PropertyTree(impl_tree.second)));

  }

  return true;

}



std::string kel::PropertyTree::getValue() const { return properties_impl_ptr_->getData<std::string>(); }
