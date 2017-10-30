//
// Created by kellerberrin on 29/10/17.
//

#ifndef KGL_VARIANT_EVIDENCE_H
#define KGL_VARIANT_EVIDENCE_H


#include "kgl_variant_db.h"
#include "kgl_genome_db.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


class VariantAnalysis {

public:

  explicit VariantAnalysis() = default;
  ~VariantAnalysis() = default;


  std::shared_ptr<GenomeVariant> SNPVariants(std::shared_ptr<ContigCountData>& count_data,
                                             std::shared_ptr<GenomeDatabase>& genome_db);

  std::shared_ptr<GenomeVariant> codonDelete(std::shared_ptr<GenomeVariant> delete_SNPs,
                                             std::shared_ptr<ContigCountData>& count_data,
                                             std::shared_ptr<GenomeDatabase>& genome_db);
private:


};



}   // namespace genome
}   // namespace kellerberrin







#endif //KGL_VARIANT_EVIDENCE_H
