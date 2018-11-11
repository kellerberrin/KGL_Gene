//
// Created by kellerberrin on 11/11/18.
//


#ifndef KGL_PROPERTIES_H
#define KGL_PROPERTIES_H

#include "kgl_exec_env.h"

#include <memory>
#include <string>


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace





// This object is a facade over the boost:: property tree object
// Properties may stored as XML or JSON. The actual format is determined at runtime
// by a file extension of "xml" or "json". Any other file extension is assumed to be
// formatted as XML. The choice of which of the two formats to use is at the
// discretion of the program user and is transparent to this object.
// The boost:: functionality is hidden using the PIMPL idiom.


class PropertyTree {

public:

  PropertyTree();
  ~PropertyTree();

  bool readProperties(const std::string& properties_file);

  bool getProperty(const std::string& property_name, std::string& property) const;

  bool getProperty(const std::string& property_name, size_t& property) const;

  // Checks for the existence of a property.
  bool checkProperty(const std::string& property_name) const;

  void treeTraversal() const;

private:

  class PropertyImpl;       // Forward declaration of the boost:: properties implementation class
  std::unique_ptr<PropertyImpl> properties_impl_ptr_;    // boost:: properties PIMPL

};




}   // namespace genome
}   // namespace kellerberrin





#endif //KGL_PROPERTIES_H
