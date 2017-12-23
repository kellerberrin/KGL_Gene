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


class VariantFactory { // Trivial top level object.

public:

  explicit VariantFactory() = default;
  virtual ~VariantFactory() = default;

  std::shared_ptr<const GenomeVariant> createVariants(const std::string &genome_name,
                                                      const std::shared_ptr<const ContigCountData> &count_data,
                                                      const std::shared_ptr<const GenomeDatabase> &genome_db,
                                                      NucleotideReadCount_t minimum_read_count,
                                                      double minimum_proportion,
                                                      Phred_t read_quality) const;

private:


};


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

  void addSingleVariant(std::shared_ptr<GenomeVariant> genome_single_variants, const Variant &variant) const;

};




}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_FACTORY_H
