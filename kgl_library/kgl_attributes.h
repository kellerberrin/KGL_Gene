//
// Created by kellerberrin on 9/11/17.
//

#ifndef KGL_ATTRIBUTES_H
#define KGL_ATTRIBUTES_H

#include <string>
#include <vector>
#include <map>


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Attributes Object to hold { key=value } pairs.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using AttributeMap = std::multimap<const std::string, std::string>;
class Attributes {

public:

  explicit Attributes() = default;
  Attributes(const Attributes&) = default;
  ~Attributes() = default;

  Attributes& operator=(const Attributes&) = default;

  // General access routines.
  bool getAttributes(const std::string &key, std::vector<std::string> &value_vec) const; // false if no key.
  void insertAttribute(const std::string& key, const std::string& value); // Always succeeds; keys are uppercase.
  void getAllAttributes(std::vector<std::pair<std::string, std::string>>& all_key_value_pairs) const;

  // Attribute keys.
  constexpr static const char* ID_KEY = "ID";
  constexpr static const char* SUPER_FEATURE_KEY = "PARENT";
  constexpr static const char* ASSIGNED_FEATURE_KEY = "ASSIGNEDFEAT";
  constexpr static const char* DESCRIPTION_KEY = "DESCRIPTION";

  // Convenience access routines.
  bool getIds(std::vector<std::string> &value_vec) const { return getAttributes(ID_KEY, value_vec); }
  bool getSuperFeatureIds(std::vector<std::string> &value_vec) const { return getAttributes(SUPER_FEATURE_KEY, value_vec); }
  bool getAssignedFeatureIds(std::vector<std::string> &value_vec) const { return getAttributes(ASSIGNED_FEATURE_KEY, value_vec); }
  bool getDescription(std::vector<std::string> &value_vec) const { return getAttributes(DESCRIPTION_KEY, value_vec); }


private:

  AttributeMap attributes_;

};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_ATTRIBUTES_H
