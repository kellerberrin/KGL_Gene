//
// Created by kellerberrin on 22/11/17.
//

#ifndef KGL_VARIANT_FACTORY_H
#define KGL_VARIANT_FACTORY_H



#include "kgl_variant_db.h"
#include "kgl_variant_compound.h"
#include "kgl_genome_db.h"
#include "kgl_variant.h"
#include "kgl_variant_factory_vcf_impl.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The top level variant factory object
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VariantFactory {

public:

  explicit VariantFactory() = default;
  virtual ~VariantFactory() = default;

  void readCountVariants(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                         std::shared_ptr<PhasedPopulation> pop_variant_ptr,
                         const std::string& genome_name,
                         const std::string& variant_file_name,
                         Phred_t read_quality,
                         NucleotideReadCount_t min_read_count,
                         double min_proportion) const;


  void readVCFVariants(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                       std::shared_ptr<UnphasedPopulation> vcf_population_ptr,
                       const std::string& genome_name,
                       const std::string& variant_file_name,
                       Phred_t variant_quality,
                       NucleotideReadCount_t min_read_count,
                       double min_proportion) const;

  size_t addGenomeSingleThreadVariant(std::shared_ptr<GenomeVariant> genome_variants, std::shared_ptr<const Variant> variant_ptr) const;

  static std::shared_ptr<const GenomeVariant> aggregateVariants(const std::shared_ptr<const GenomeDatabase>& genome_db_ptr,
                                                                const std::string& genome_name,
                                                                const std::shared_ptr<const GenomeVariant>& single_variant_ptr);

private:

  bool isFileNamePrefix(const std::string& prefix, const std::string& variant_file_name) const;

  void addGenome(std::shared_ptr<const GenomeVariant> genome_variant_ptr,
                 std::shared_ptr<PhasedPopulation> pop_variant_ptr,
                 Phred_t read_quality) const;

  std::shared_ptr<const GenomeVariant> createSamVariants(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                                         const std::string& genome_name,
                                                         const std::string& sam_file_name,
                                                         Phred_t read_quality,
                                                         NucleotideReadCount_t min_read_count,
                                                         double min_proportion) const;

  std::shared_ptr<const GenomeVariant> createBamVariants(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                                         const std::string& genome_name,
                                                         const std::string& sam_file_name,
                                                         Phred_t read_quality,
                                                         NucleotideReadCount_t min_read_count,
                                                         double min_proportion) const;


  constexpr static const char* SAM_FILE_EXTENSTION_ = ".SAM";
  constexpr static const char* BAM_FILE_EXTENSTION_ = ".BAM";
  constexpr static const char* VCF_FILE_EXTENSTION_ = ".VCF";
  constexpr static const char* GATK_FILE_PREFIX_ = "GATK";
  constexpr static const char* PF3K_FILE_PREFIX_ = "PF3K";


};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_FACTORY_H
