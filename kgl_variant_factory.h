//
// Created by kellerberrin on 22/11/17.
//

#ifndef KGL_VARIANT_FACTORY_H
#define KGL_VARIANT_FACTORY_H



#include "kgl_variant_db.h"
#include "kgl_variant_compound.h"
#include "kgl_genome_db.h"
#include "kgl_variant.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The top level variant factory object
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VariantFactory {

public:

  explicit VariantFactory() = default;
  virtual ~VariantFactory() = default;

  std::shared_ptr<const GenomeVariant> createVariants(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                                      const std::string& genome_name,
                                                      const std::string& variant_file_name,
                                                      Phred_t read_quality,
                                                      Phred_t variant_quality,
                                                      NucleotideReadCount_t min_read_count,
                                                      double min_proportion) const;

  // Annotates variants with genome information and adds to the variant genome. For use in the SAM, BAM and VCF parsers.
  static size_t addSingleVariant(std::shared_ptr<GenomeVariant> genome_single_variants, const Variant &variant);

private:

  std::shared_ptr<const GenomeVariant> createSamVariants(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                                         const std::string& genome_name,
                                                         const std::string& sam_file_name,
                                                         Phred_t read_quality,
                                                         Phred_t variant_quality,
                                                         NucleotideReadCount_t min_read_count,
                                                         double min_proportion) const;

  std::shared_ptr<const GenomeVariant> createBamVariants(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                                         const std::string& genome_name,
                                                         const std::string& sam_file_name,
                                                         Phred_t read_quality,
                                                         Phred_t variant_quality,
                                                         NucleotideReadCount_t min_read_count,
                                                         double min_proportion) const;

  std::shared_ptr<const GenomeVariant> createVcfVariants(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                                         const std::string& genome_name,
                                                         const std::string& vcf_file_name,
                                                         Phred_t read_quality,
                                                         Phred_t variant_quality,
                                                         NucleotideReadCount_t min_read_count,
                                                         double min_proportion) const;

  std::shared_ptr<const GenomeVariant> aggregateVariants(const std::shared_ptr<const GenomeDatabase>& genome_db_ptr,
                                                         const std::string& genome_name,
                                                         const std::shared_ptr<const GenomeVariant>& single_variant_ptr) const;

  constexpr static const char* SAM_FILE_EXTENSTION_ = ".SAM";
  constexpr static const char* BAM_FILE_EXTENSTION_ = ".BAM";
  constexpr static const char* VCF_FILE_EXTENSTION_ = ".VCF";

};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_FACTORY_H
