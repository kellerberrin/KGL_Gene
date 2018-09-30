//
// Created by kellerberrin on 19/03/18.
//

#ifndef KGL_GENOME_AUX_CSV_H
#define KGL_GENOME_AUX_CSV_H


#include <string>
#include <vector>
#include <map>
#include "kgl_genome_types.h"
#include "kgl_variant_db.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Reads up the Pf3k auxiliary data description file in csv format, 1 line for each sample.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using CountryPair = std::pair<GenomeId_t, std::shared_ptr<const PhasedPopulation>>;

using AuxAttributeVector = std::vector<std::string>;
using AuxSampleMap = std::map<std::string, AuxAttributeVector>;


class GenomeAuxData {

public:

  GenomeAuxData() = default;
  ~GenomeAuxData() = default;

// Convenience function splits the phased populations into different countries (preferred genomes only).
  static std::vector<CountryPair> getCountries(const std::string& aux_file, std::shared_ptr<const PhasedPopulation> population_ptr);

  bool readParseAuxData(const std::string& aux_file_name);

  std::string locationDate(const std::string& genome_name) const;

  bool isFieldSample(const GenomeId_t& genome_name) const;

  bool isPreferredSample(const GenomeId_t& genome_name) const;

  bool locationDate(const std::string& genome_name,
                    std::string& country,
                    std::string& location,
                    std::string& year) const;

  const AuxSampleMap& getMap() const { return aux_sample_information_; }

  std::vector<std::string> countryList() const;
  std::vector<GenomeId_t> getCountry(const std::string& country) const;
  std::vector<GenomeId_t> getFieldSamples() const;
  std::vector<GenomeId_t> getPreferredSamples() const;

private:

  AuxAttributeVector aux_data_header_;  // Always assumed to be the first line (uppercase, query case conversion automatic).
  AuxSampleMap aux_sample_information_;

  bool parseHeader(const std::string& record_str);
  bool parseDataline(const std::string& record_str);
  bool tokenize(const std::string& parse_text, AuxAttributeVector& attribute_vector);
  size_t fieldOffset(const std::string& field_name) const;
  bool getGenomeAttributes(const std::string& genome_name, AuxAttributeVector& attribute_vector) const;

  constexpr static const char* FIELD_SAMPLE = "IsFieldSample";
  constexpr static const char* PREFERRED_SAMPLE = "PreferredSample";
  constexpr static const char* COUNTRY = "country";
  constexpr static const char* SITE = "site";
  constexpr static const char* COLLECTION_YEAR = "collection_year";
  constexpr static const char* TRUE_FLAG = "TRUE";

};






}   // namespace genome
}   // namespace kellerberrin




#endif //KGL_KGL_GENOME_AUX_CSV_H
