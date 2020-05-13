//
// Created by kellerberrin on 26/12/17.
//

#ifndef KGL_VARIANT_FACTORY_VCF_H
#define KGL_VARIANT_FACTORY_VCF_H

#include <memory>
#include <string>
#include <map>
#include "kgl_variant_db.h"


namespace kellerberrin::genome {   //  organization level namespace


class VcfFactory  {

public:

  VcfFactory() = default;
  ~VcfFactory() = default;

  [[nodiscard]] static std::shared_ptr<UnphasedPopulation>
  gatkMultiGenomeVCFVariants( const std::shared_ptr<const GenomeReference> genome_db_ptr,
                              const std::string &vcf_file_name,
                              const EvidenceInfoSet& evidence_set);

  [[nodiscard]] static std::shared_ptr<UnphasedGenome>
  GRChNoGenomeVCFVariants( const std::shared_ptr<const GenomeReference> genome_db_ptr,
                           const std::string &vcf_file_name,
                           const ContigAliasMap& contig_alias_map,
                           const EvidenceInfoSet& evidence_set);


private:



};


}   // end namespace


#endif //KGL_VARIANT_FACTORY_VCF_H
