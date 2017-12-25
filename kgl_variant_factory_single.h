//
// Created by kellerberrin on 25/12/17.
//

#ifndef KGL_VARIANT_FACTORY_SINGLE_H
#define KGL_VARIANT_FACTORY_SINGLE_H



#include "kgl_variant_db.h"
#include "kgl_variant_compound.h"
#include "kgl_genome_db.h"
#include "kgl_variant.h"
#include "kgl_variant_factory.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Creates read count based SNP variants.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class SingleFactory : public VariantFactory {

public:

  explicit SingleFactory() = default;
  ~SingleFactory() override = default;


  std::shared_ptr<const GenomeVariant> createSingleVariants(const std::string &genome_name,
                                                            const std::shared_ptr<const ContigCountData> &count_data,
                                                            const std::shared_ptr<const GenomeDatabase> &genome_db_ptr,
                                                            NucleotideReadCount_t minimum_read_count,
                                                            double minimum_proportion,
                                                            Phred_t read_quality) const;

private:

  size_t GenerateSNPDelete(const std::string &genome_name,
                           std::shared_ptr<ContigFeatures> contig_ptr,
                           ContigOffset_t contig_offset,
                           DNA5::Alphabet reference_nucleotide,
                           const NucleotideReadCount_t nucleotide_count_array[],
                           NucleotideReadCount_t minimum_read_count,
                           double minimum_proportion,
                           Phred_t read_quality,
                           std::shared_ptr<GenomeVariant> genome_single_variants) const;

  size_t GenerateInsert(const std::string &genome_name,
                        std::shared_ptr<ContigFeatures> contig_ptr,
                        ContigOffset_t contig_offset,
                        DNA5::Alphabet reference_nucleotide,
                        const NucleotideReadCount_t nucleotide_count_array[],
                        NucleotideReadCount_t minimum_read_count,
                        double minimum_proportion,
                        Phred_t read_quality,
                        std::shared_ptr<GenomeVariant> genome_single_variants) const;

  size_t addSingleVariant(std::shared_ptr<GenomeVariant> genome_single_variants, const Variant &variant) const;

};




}   // namespace genome
}   // namespace kellerberrin



#endif //READSAMFILE_KGL_VARIANT_FACTORY_SINGLE_H
