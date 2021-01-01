//
// Created by kellerberrin on 15/2/20.
//

#ifndef KGL_PROPERTY_TREE_H
#define KGL_PROPERTY_TREE_H


#include "kel_exec_env.h"

#include <memory>
#include <string>
#include <set>


namespace kellerberrin {  //  organization level namespace


// This object is a facade over the boost:: property tree object
// Properties may stored as XML or JSON. The actual format is determined at runtime
// by a file extension of "xml" or "json". Any other file extension is assumed to be
// formatted as XML. The choice of which of the two formats to use is at the
// discretion of the program user and is transparent to this object.
// The boost:: functionality is hidden using the PIMPL idiom.

class PropertyTree; // fwd
using SubPropertyTree = std::pair<std::string, PropertyTree>;

class PropertyTree {

public:

  class PropertyImpl;       // Forward declaration of the boost:: properties implementation class

  PropertyTree(); // defined in implementation file
  PropertyTree(const PropertyTree&);
  PropertyTree(const PropertyImpl&);

  ~PropertyTree(); // defined in implementation file

  [[nodiscard]] bool readProperties(const std::string& properties_file);

  [[nodiscard]] bool getProperty(const std::string& property_name, std::string& property) const;

  [[nodiscard]] bool getOptionalProperty(const std::string& property_name, std::string& property) const;

  [[nodiscard]] bool getPropertyVector(const std::string& property_name, std::vector<std::string>& property_vector) const;

  [[nodiscard]] bool getNodeVector(const std::string& node_name, std::vector<std::string>& node_vector) const;

  [[nodiscard]] bool getPropertyTreeVector(const std::string& property_name, std::vector<SubPropertyTree>& property_tree_vector) const;

  [[nodiscard]] bool getPropertySubTreeVector(std::vector<SubPropertyTree>& property_tree_vector) const;

  [[nodiscard]] bool getProperty(const std::string& property_name, size_t& property) const;

  [[nodiscard]] bool getFileProperty(const std::string& property_name, const std::string& work_directory, std::string& file_path) const;

  [[nodiscard]] bool getOptionalFileProperty(const std::string& property_name, const std::string& work_directory, std::string& file_path) const;

  // Checks for the existence of a property.
  [[nodiscard]] bool checkProperty(const std::string& property_name) const;

  [[nodiscard]] std::string getValue() const;

  void treeTraversal() const;

private:

  std::unique_ptr<PropertyImpl> properties_impl_ptr_;    // boost:: properties PIMPL


};


}   // end namespace




#endif //KGL_KGL_PROPERTY_TREE_H
