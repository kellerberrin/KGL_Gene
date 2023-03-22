//
// Created by kellerberrin on 9/11/17.
//

#ifndef KGL_ATTRIBUTES_H
#define KGL_ATTRIBUTES_H

#include <string>
#include <vector>
#include <map>


namespace kellerberrin::genome {   //  organization::project level namespace


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
  [[nodiscard]] std::vector<std::string> getAttributes(const std::string &key) const; // false if no key.
  void insertAttribute(const std::string& key, const std::string& value); // Always succeeds; keys are uppercase.
  void insertAttribute(std::string&& key, std::string&& value);

  // Attribute keys.
  constexpr static const char* ID_KEY = "ID";
  constexpr static const char* SUPER_FEATURE_KEY = "PARENT";
  constexpr static const char* ASSIGNED_FEATURE_KEY = "ASSIGNEDFEAT";
  constexpr static const char* DESCRIPTION_KEY = "DESCRIPTION";
  constexpr static const char* NAME_KEY = "NAME";
  constexpr static const char* GENE_BIOTYPE_KEY = "GENE_BIOTYPE";
  constexpr static const char SUPER_FEATURE_DELIMITER = ',';

  // Convenience access routines.
  [[nodiscard]] bool getIds(std::vector<std::string> &value_vec) const { return getAttributes(ID_KEY, value_vec); }
  // Super features are further parsed by comma ','. For example; "PF3D7_0108400.1,PF3D7_0108400.2" returns 2 super features.
  void getSuperFeatureIds(std::vector<std::string> &value_vec) const;
  void getAssignedFeatureIds(std::vector<std::string> &value_vec) const { getAttributes(ASSIGNED_FEATURE_KEY, value_vec); }
  bool getDescription(std::vector<std::string> &value_vec) const { return getAttributes(DESCRIPTION_KEY, value_vec); }
  [[nodiscard]] std::vector<std::string> getDescription() const { return getAttributes(DESCRIPTION_KEY); }
  bool getGeneBioType(std::vector<std::string> &value_vec) const { return getAttributes(GENE_BIOTYPE_KEY, value_vec); }
  bool getName(std::vector<std::string> &value_vec) const { return getAttributes(NAME_KEY, value_vec); }
  [[nodiscard]] const AttributeMap& getMap() const { return attributes_; }

  // Get HGNC identifiers from Homo Sapien GFF files.
  std::string getHGNC() const;

  // Mainly used for tesing.
  [[nodiscard]] bool equivalent(const Attributes& lhs) const { return attributes_ == lhs.attributes_; }

private:

  AttributeMap attributes_;

  constexpr static const char* DBXREF_ = "DBXREF";
  constexpr static const char* HGNC_ = "HGNC:";

};


}   // end namespace


#endif //KGL_ATTRIBUTES_H
