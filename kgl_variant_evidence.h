//
// Created by kellerberrin on 29/10/17.
//

#ifndef KGL_VARIANT_EVIDENCE_H
#define KGL_VARIANT_EVIDENCE_H


#include "kgl_variant_db.h"
#include "kgl_variant_compound.h"
#include "kgl_genome_db.h"
#include "kgl_variant.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


class VariantAnalysis {

public:

  explicit VariantAnalysis() = default;
  ~VariantAnalysis() = default;


  std::shared_ptr<const GenomeVariant> SNPVariants(const std::shared_ptr<const ContigCountData>& count_data,
                                                   const std::shared_ptr<const GenomeDatabase>& genome_db,
                                                   NucleotideReadCount_t minimum_read_count,
                                                   double minimum_proportion);

  std::shared_ptr<const GenomeVariant> codonDelete(const std::shared_ptr<const GenomeVariant>& delete_SNPs,
                                                   const std::shared_ptr<const ContigCountData>& count_data,
                                                   const std::shared_ptr<const GenomeDatabase>& genome_db);

  std::shared_ptr<const GenomeVariant> compoundSNP(const std::shared_ptr<const GenomeVariant>& SNPs,
                                                   const std::shared_ptr<const GenomeDatabase>& genome_db);

private:


  void aggregateCodingDeletions(const std::shared_ptr<const GenomeVariant>& delete_SNPs,
                                const std::shared_ptr<const GenomeDatabase>& genome_db_ptr,
                                std::vector<CompoundVariantMap>& contiguous_delete_vec);

  void generateCodonDeletes(const std::shared_ptr<const GenomeDatabase>& genome_db_ptr,
                            const std::shared_ptr<const ContigCountData>& count_data,
                            const std::vector<CompoundVariantMap>& contiguous_delete_vec,
                            std::shared_ptr<GenomeVariant> genome_variant_ptr);

  std::shared_ptr<const Variant> createCompoundDelete(const CompoundVariantMap& variant_map);

};



}   // namespace genome
}   // namespace kellerberrin







#endif //KGL_VARIANT_EVIDENCE_H
