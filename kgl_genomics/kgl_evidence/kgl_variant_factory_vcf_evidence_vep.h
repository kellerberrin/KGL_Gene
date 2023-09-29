//
// Created by kellerberrin on 22/04/23.
//

#ifndef KGL_VARIANT_FACTORY_VCF_EVIDENCE_VEP_H
#define KGL_VARIANT_FACTORY_VCF_EVIDENCE_VEP_H


#include "kgl_variant_factory_vcf_evidence_analysis.h"


namespace kellerberrin::genome {   //  organization level namespace




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using VepFieldValueMap = std::map<std::string, size_t>;

class VepSubFieldValues {

public:

  VepSubFieldValues(std::string vep_sub_field) : vep_sub_field_(std::move(vep_sub_field)) {}
  ~VepSubFieldValues() = default;

  [[nodiscard]] const VepFieldValueMap& getMap() const { return field_value_map_; }

  // For a single variant
  bool getSubFieldValues(const std::shared_ptr<const Variant>& variant_ptr);
  // For a population.
  bool getPopulationValues(const std::shared_ptr<const PopulationDB>& population_ptr);
  // For a genome
  bool getGenomeValues(const std::shared_ptr<const GenomeDB>& genome_ptr);
  // For a contig_ref_ptr.
  bool getContigValues(const std::shared_ptr<const ContigDB>& contig_ptr);

private:

  const std::string vep_sub_field_;
  VepFieldValueMap field_value_map_;


}; // struct.




} // namespace

#endif //KGL_VARIANT_FACTORY_VCF_EVIDENCE_VEP_H
