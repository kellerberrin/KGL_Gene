//
// Created by kellerberrin on 11/11/18.
//

#include "kgl_utility.h"
#include "kgl_properties.h"
#include "kgl_table_impl.h"


#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/tokenizer.hpp>


namespace kgl = kellerberrin::genome;
namespace pt = boost::property_tree;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PropertyTree::PropertyImpl does all the heavy lifting using boost::property_tree.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class kgl::PropertyTree::PropertyImpl {

public:

  PropertyImpl() = default;
  ~PropertyImpl() = default;

  bool readPropertiesFile(const std::string& properties_file);

  bool getProperty(const std::string& property_name, std::string& property) const;

  bool getPropertyVector(const std::string& property_name, std::vector<std::string>& property_vector) const;

  bool getProperty(const std::string& property_name, size_t& property) const;

  void treeTraversal() const;

  bool checkProperty(const std::string& property_name) const;


private:

  constexpr static const char* JSON_FILE_EXT_ = "JSON";

  // boost property tree object
  pt::ptree property_tree_;

  void printTree(const std::string& parent, const pt::ptree& property_tree) const;

};



bool kgl::PropertyTree::PropertyImpl::readPropertiesFile(const std::string& properties_file) {

  std::string uc_file_extenstion = Utility::toupper(Utility::fileExtension(properties_file));

  try {

    if (uc_file_extenstion == JSON_FILE_EXT_) {

      pt::read_json(properties_file, property_tree_);

    } else {

      pt::read_xml(properties_file, property_tree_);

    }

  }
  catch(...) {

    ExecEnv::log().error("PropertyTree; Missing or Malformed property tree in file: {}", properties_file);
    return false;

  }

  return true;

}


bool kgl::PropertyTree::PropertyImpl::checkProperty(const std::string& property_name) const {

  try {

    property_tree_.get<std::string>(property_name);

  }
  catch (...) {

    return false;

  }

  return true;

}


bool kgl::PropertyTree::PropertyImpl::getProperty(const std::string& property_name, std::string& property) const {

  try {

    property = property_tree_.get<std::string>(property_name);
    property = Utility::trimEndWhiteSpace(property);

  }
  catch (...) {

    ExecEnv::log().error("PropertyTree; Property: {} not found", property_name);
    ExecEnv::log().error("***********Property Tree Contents*************");
    treeTraversal();
    ExecEnv::log().error("**********************************************");
    return false;

  }

  return true;

}


bool kgl::PropertyTree::PropertyImpl::getPropertyVector(const std::string& property_name, std::vector<std::string>& property_vector) const {

  try {

    for (auto child : property_tree_.get_child(property_name)) {
      // The data function is used to access the data stored in a node.
      property_vector.push_back(Utility::trimEndWhiteSpace(child.second.data()));

    }

  }
  catch (...) {

    ExecEnv::log().error("PropertyTree; Property: {} not found", property_name);
    ExecEnv::log().error("***********Property Tree Contents*************");
    treeTraversal();
    ExecEnv::log().error("**********************************************");
    return false;

  }

  return true;

}



bool kgl::PropertyTree::PropertyImpl::getProperty(const std::string& property_name, size_t& property) const {

  std::string size_string;
  if (not getProperty(property_name, size_string)) {

    return false;

  }

  try {

    property= std::stoull(size_string);

  }
  catch (...) {

    ExecEnv::log().error("PropertyTree; Property: {}, Value: {} is not an unsigned integer", property_name, size_string);
    return false;

  }

  return true;

}


void kgl::PropertyTree::PropertyImpl::treeTraversal() const {

  printTree("", property_tree_);

}

void kgl::PropertyTree::PropertyImpl::printTree(const std::string& parent, const pt::ptree& property_tree) const {

  for (auto item : property_tree) {

    std::string key = item.first;
    std::string parent_key;
    if (parent.empty()) {

      parent_key = key;

    } else {

      parent_key = parent + "." + key;

    }
    std::string value = property_tree.get<std::string>(key);
    value = Utility::trimEndWhiteSpace(value);
    ExecEnv::log().info("PropertyTree; key: {}, value: {}", parent_key, value);
    printTree(parent_key, item.second);

  }

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PropertyTree is a public facade class that passes the functionality onto PropertyTree::PropertyImpl.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::PropertyTree::PropertyTree() : properties_impl_ptr_(std::make_unique<kgl::PropertyTree::PropertyImpl>()) {}
kgl::PropertyTree::~PropertyTree() {}  // DO NOT DELETE or USE DEFAULT. Required because of incomplete pimpl type.

// Functionality passed to the implmentation.

bool kgl::PropertyTree::readProperties(const std::string& properties_file) {

  return properties_impl_ptr_->readPropertiesFile(properties_file);

}


bool kgl::PropertyTree::getProperty(const std::string& property_name, std::string& property) const {

  return properties_impl_ptr_->getProperty(property_name, property);

}


bool kgl::PropertyTree::getProperty(const std::string& property_name, size_t& property) const {

  return properties_impl_ptr_->getProperty(property_name, property);

}

void kgl::PropertyTree::treeTraversal() const {

  return properties_impl_ptr_->treeTraversal();

}

bool kgl::PropertyTree::checkProperty(const std::string& property_name) const {

  return properties_impl_ptr_->checkProperty(property_name);

}


bool kgl::PropertyTree::getFileProperty(const std::string& property_name, const std::string& work_directory, std::string& file_path) const {

  std::string file_name;

  if (not getProperty(property_name, file_name)) {

    ExecEnv::log().warn("PropertyTree::getFileProperty; Requested file property: {} not found. A list of all valid properties follows:",
                        property_name);

    treeTraversal();
    return false;

  }

  file_path = Utility::filePath(file_name, work_directory);

  return Utility::fileExists(file_path);

}


bool kgl::PropertyTree::getPropertyVector(const std::string& property_name, std::vector<std::string>& property_vector) const {

  return properties_impl_ptr_->getPropertyVector(property_name, property_vector);

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// High level application specific property retrieval
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::RuntimeProperties::readProperties(const std::string& properties_file) {

  std::string properties_path = Utility::filePath(properties_file, work_directory_);
  return property_tree_.readProperties(properties_path);

}



void kgl::RuntimeProperties::getGenomeDBFiles(const std::string& organism,
                                              std::string& fasta_file,
                                              std::string& gff_file,
                                              std::string& gaf_file,
                                              std::string& translationTable) const {


  const std::string dot = ".";

  std::string key = std::string(RUNTIME_ROOT_) + organism + dot + std::string(GAF_FILE_) + std::string(VALUE_);
  property_tree_.getFileProperty(key, work_directory_, gaf_file);

  key = std::string(RUNTIME_ROOT_) + organism + dot + std::string(FASTA_FILE_) + std::string(VALUE_);
  if (not property_tree_.getFileProperty(key, work_directory_, fasta_file)) {

    ExecEnv::log().critical("Required Fasta File: {} does not exist.", fasta_file);

  }

  key = std::string(RUNTIME_ROOT_) + organism + dot + std::string(GFF_FILE_) + std::string(VALUE_);
  if (not property_tree_.getFileProperty(key, work_directory_, gff_file)) {

    ExecEnv::log().critical("Required GFF3 File: {} does not exist.", gff_file);

  }

  key = std::string(RUNTIME_ROOT_) + organism + dot + std::string(TRANSLATION_TABLE_) + std::string(VALUE_);
  if (not property_tree_.getProperty(key, translationTable)) {

    translationTable = Tables::STANDARDTABLE->table_name; // The standard amino acid translation table.
    ExecEnv::log().warn("Amino Acid translation table not specified, using standard table: {}", translationTable);

  }

}



bool kgl::RuntimeProperties::getGenomeAuxFiles(const std::string& organism,
                                               const std::string& aux_file_name,
                                               std::string& aux_file) const {

  const std::string dot = ".";

  std::string key = std::string(RUNTIME_ROOT_)  + organism + dot + aux_file_name + std::string(VALUE_);

  return property_tree_.getFileProperty(key, work_directory_, aux_file);

}


size_t kgl::RuntimeProperties::getVCFPloidy() const {

  std::string key = std::string(RUNTIME_ROOT_) + std::string(VCF_PLOIDY_) + std::string(VALUE_);
  size_t ploidy;

  if (not property_tree_.getProperty(key, ploidy)) {

    ExecEnv::log().warn("VCF ploidy is not defined (need not be organism ploidy)");
    ploidy = DEFAULT_PLOIDY_;

  }

  return ploidy;

}


bool kgl::RuntimeProperties::getMixtureFile(std::string& mixture_file) const {

  std::string key = std::string(RUNTIME_ROOT_) + std::string(MIXTURE_FILE_) + std::string(VALUE_);
  if (not property_tree_.getProperty(key, mixture_file)) {

    return false;

  }

  mixture_file = Utility::filePath(mixture_file, work_directory_);

  return Utility::fileExists(mixture_file);

}


bool kgl::RuntimeProperties::getAuxFile(std::string& aux_file) const {

  std::string key = std::string(RUNTIME_ROOT_) + std::string(AUX_FILE_) + std::string(VALUE_);
  if (not property_tree_.getProperty(key, aux_file)) {

    return false;

  }

  aux_file = Utility::filePath(aux_file, work_directory_);

  return Utility::fileExists(aux_file);

}


bool kgl::RuntimeProperties::getVCFFiles(std::vector<std::string>& vcf_files) const {

  std::string key = std::string(RUNTIME_ROOT_) + std::string(VCF_LIST_) + std::string(FILE_LIST_);
  std::vector<std::string> truncated_vcf_vector;
  if (not property_tree_.getPropertyVector(key, truncated_vcf_vector)) {

    return false;

  }

  for (auto file : truncated_vcf_vector) {

    vcf_files.push_back(Utility::filePath(file, work_directory_));

  }

  return true;

}

