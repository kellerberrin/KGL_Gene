//
// Created by kellerberrin on 29/10/17.
//

#ifndef KGL_VARIANT_EVIDENCE_H
#define KGL_VARIANT_EVIDENCE_H


#include "kgl_variant_db.h"
#include "kgl_variant_compound.h"
#include "kgl_genome_db.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


class VariantAnalysis {

public:

  explicit VariantAnalysis() = default;
  ~VariantAnalysis() = default;

  std::shared_ptr<const GenomeVariant> SNPVariants(std::shared_ptr<const ContigCountData> count_data,
                                                   std::shared_ptr<const GenomeDatabase> genome_db);

  std::shared_ptr<const GenomeVariant> codonDelete(std::shared_ptr<const GenomeVariant> delete_SNPs,
                                                   std::shared_ptr<const ContigCountData> count_data,
                                                   std::shared_ptr<const GenomeDatabase> genome_db);
private:


  void aggregateCodingDeletions(std::shared_ptr<const GenomeVariant>& delete_SNPs,
                                std::shared_ptr<const GenomeDatabase>& genome_db_ptr,
                                std::vector<CompoundVariantMap>& contiguous_delete_vec);

  bool membershipCodingDeletions(const std::vector<CompoundVariantMap>& contiguous_delete_vec);

  void generateCodonDeletes(std::shared_ptr<const GenomeDatabase>& genome_db_ptr,
                            std::shared_ptr<const ContigCountData> count_data,
                            const std::vector<CompoundVariantMap>& contiguous_delete_vec);
  void printCompoundVariant(const std::vector<CompoundVariantMap>& contiguous_delete_vec);
};



}   // namespace genome
}   // namespace kellerberrin







#endif //KGL_VARIANT_EVIDENCE_H
