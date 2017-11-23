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
  ~VariantFactory() = default;


private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Creates read count based SNP variants.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VariantSNPFactory : public VariantFactory {

public:

  explicit VariantSNPFactory() = default;
  ~VariantSNPFactory() = default;


  std::shared_ptr<const GenomeVariant> SNPVariants(const std::string& genome_name,
                                                   const std::shared_ptr<const ContigCountData>& count_data,
                                                   const std::shared_ptr<const GenomeDatabase>& genome_db,
                                                   NucleotideReadCount_t minimum_read_count,
                                                   double minimum_proportion);

private:

  void addSNPVariant(std::shared_ptr<GenomeVariant> genome_snp_variants, const SNPVariantDNA5& variant);

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The top level compound variant factory.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VariantCompoundFactory : public VariantFactory {

public:

  explicit VariantCompoundFactory() = default;
  ~VariantCompoundFactory() = default;


  std::shared_ptr<GenomeVariant>
  disaggregateCompoundVariants(const std::shared_ptr<const GenomeVariant>& genome_variant,
                               const std::shared_ptr<const GenomeDatabase>& genome_db) const;



private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The compound delete variant factory.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VariantDeleteFactory : public VariantCompoundFactory {

public:

  explicit VariantDeleteFactory() = default;
  ~VariantDeleteFactory() = default;


  std::shared_ptr<const GenomeVariant> codonDelete(const std::shared_ptr<const GenomeVariant>& delete_SNPs,
                                                   const std::shared_ptr<const ContigCountData>& count_data,
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



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The compound SNP variant factory.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VariantCompoundSNPFactory : public VariantCompoundFactory {

public:

  explicit VariantCompoundSNPFactory() = default;
  ~VariantCompoundSNPFactory() = default;


  std::shared_ptr<const GenomeVariant> compoundSNP(const std::shared_ptr<const GenomeVariant>& SNPs,
                                                   const std::shared_ptr<const GenomeDatabase>& genome_db);

private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The compound insert variant factory.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VariantInsertFactory : public VariantCompoundFactory {

public:

  explicit VariantInsertFactory() = default;
  ~VariantInsertFactory() = default;



private:


};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_FACTORY_H
