//
// Created by kellerberrin on 11/11/18.
//

#include "kgl_utility.h"
#include "kgl_properties.h"


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

    ExecEnv::log().error("PropertyTree; Malformed property tree in file: {}", properties_file);
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
