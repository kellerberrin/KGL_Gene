//
// Created by kellerberrin on 15/2/20.
//

#ifndef KGL_PROPERTY_TREE_H
#define KGL_PROPERTY_TREE_H


#include "kel_exec_env.h"

#include <memory>
#include <string>
#include <shared_mutex>
#include <utility>
#include <vector>


namespace kellerberrin {  //  organization level namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This object is a facade over the boost:: property tree object
// Properties may stored as XML or JSON. The actual format is determined at runtime
// by a file extension of "xml" or "json". Any other file extension is assumed to be
// formatted as XML. The choice of which of the two formats to use is at the
// discretion of the program user and is transparent to this object.
// The boost:: functionality is hidden using the PIMPL idiom.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class PropertyTree; // fwd
using SubPropertyTree = std::pair<std::string, PropertyTree>;

/// PropertyTree is a read-only facade over the boost::property_tree object, hiding the
/// underlying boost::implementation via the PIMPL idiom. Parsed trees may be shared by
/// multiple reader threads; accessor methods take a shared lock, and readProperties()
/// atomically swaps in a newly parsed tree under a unique lock.
class PropertyTree {

public:

  class PropertyImpl;       // Forward declaration of the boost:: properties implementation class.

  PropertyTree(); // Defined in implementation file.
  PropertyTree(const PropertyTree&);
  explicit PropertyTree(const PropertyImpl&);

  ~PropertyTree(); // Defined in implementation file.

  /// Reads and parses an XML (default) or JSON properties file. Returns true on success.
  [[nodiscard]] bool readProperties(const std::string& properties_file);

  /// Retrieves the string value of the named property. Returns false if not present.
  [[nodiscard]] bool getProperty(const std::string& property_name, std::string& property) const;

  /// Retrieves the optional string property. Returns false if absent (never an error).
  [[nodiscard]] bool getOptionalProperty(const std::string& property_name, std::string& property) const;

  /// Retrieves a vector of string values from a named node.
  [[nodiscard]] bool getPropertyVector(const std::string& property_name, std::vector<std::string>& property_vector) const;

  /// Retrieves a vector of node data values matching the node name.
  [[nodiscard]] bool getNodeVector(const std::string& node_name, std::vector<std::string>& node_vector) const;

  /// Retrieves a vector of sub-trees from the named node.
  [[nodiscard]] bool getPropertyTreeVector(const std::string& property_name, std::vector<SubPropertyTree>& property_tree_vector) const;

  /// Retrieves a vector of all immediate sub-trees.
  [[nodiscard]] bool getPropertySubTreeVector(std::vector<SubPropertyTree>& property_tree_vector) const;

  /// Reads a numeric (size_t) property.
  [[nodiscard]] bool getProperty(const std::string& property_name, size_t& property) const;

  /// Reads a file path property relative to work_directory and validates that the file exists.
  [[nodiscard]] bool getFileProperty(const std::string& property_name, const std::string& work_directory, std::string& file_path) const;

  /// Reads a file path property, creating a zero-length file if it does not already exist.
  [[nodiscard]] bool getFileCreateProperty(const std::string& property_name, const std::string& work_directory, std::string& file_path) const;

  /// Reads an optional file path property. Returns false if the property is absent.
  [[nodiscard]] bool getOptionalFileProperty(const std::string& property_name, const std::string& work_directory, std::string& file_path) const;

  /// Checks for the existence of a property.
  [[nodiscard]] bool checkProperty(const std::string& property_name) const;

  /// Returns the string value of the root node.
  [[nodiscard]] std::string getValue() const;

  /// Recursively logs the entire tree contents.
  void treeTraversal() const;

private:

  mutable std::shared_mutex tree_mutex_;                   // Synchronizes reads and replacement.
  std::unique_ptr<PropertyImpl> properties_impl_ptr_;      // boost:: properties PIMPL

  // Ignore these.
  constexpr static const char HELP_[] = "help";
  constexpr static const char COMMENT_[] = "<xmlcomment>";

};


}   // end namespace




#endif //KGL_KGL_PROPERTY_TREE_H
