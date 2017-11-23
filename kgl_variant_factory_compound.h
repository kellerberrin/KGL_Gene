//
// Created by kellerberrin on 23/11/17.
//

#ifndef KGL_VARIANT_FACTORY_COMPOUND_H
#define KGL_VARIANT_FACTORY_COMPOUND_H


#include "kgl_variant_db.h"
#include "kgl_variant_compound.h"
#include "kgl_genome_db.h"
#include "kgl_variant_snp.h"
#include "kgl_variant_factory.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The top level compound variant factory.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VariantCompoundFactory : public VariantFactory {

public:

  explicit VariantCompoundFactory() = default;
  ~VariantCompoundFactory() override = default;

  // Disaggregate compound variant.
  std::shared_ptr<GenomeVariant>
  disaggregateCompoundVariants(const std::shared_ptr<const GenomeVariant>& genome_variant,
                               const std::shared_ptr<const GenomeDatabase>& genome_db) const;

protected:

  // Generate compound variant by aggregating consecutive variants in a coding sequence (see implementation for detail).
  bool aggregateCompoundVariants(const std::shared_ptr<const GenomeVariant>& genome_variants,
                                 std::vector<std::shared_ptr<const CompoundVariantMap>>& aggregated_variants_vec) const;

  // Determine which variants to aggregate.
  virtual bool aggregateVariant(const std::shared_ptr<const Variant>& variant_ptr) const = 0;

private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The compound delete variant factory.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VariantDeleteFactory : public VariantCompoundFactory {

public:

  explicit VariantDeleteFactory() = default;
  ~VariantDeleteFactory() override = default;


  std::shared_ptr<const GenomeVariant> compoundDelete(const std::shared_ptr<const GenomeVariant>& variants,
                                                      const std::shared_ptr<const GenomeDatabase>& genome_db);


private:

  bool aggregateVariant(const std::shared_ptr<const Variant>& variant_ptr) const override;

  void generateCompoundDeletes( const std::vector<std::shared_ptr<const CompoundVariantMap>>& contiguous_delete_vec,
                                std::shared_ptr<GenomeVariant>& genome_variant_ptr);

  std::shared_ptr<const Variant> createCompoundDelete(const CompoundVariantMap& variant_map);

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The compound insert variant factory.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VariantInsertFactory : public VariantCompoundFactory {

public:

  explicit VariantInsertFactory() = default;
  ~VariantInsertFactory() = default;

  std::shared_ptr<const GenomeVariant> compoundInsert(const std::shared_ptr<const GenomeVariant>& SNPs,
                                                      const std::shared_ptr<const GenomeDatabase>& genome_db);

private:

  bool aggregateVariant(const std::shared_ptr<const Variant>& variant_ptr) const override;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The compound SNP variant factory.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VariantCompoundSNPFactory : public VariantCompoundFactory {

public:

  explicit VariantCompoundSNPFactory() = default;
  ~VariantCompoundSNPFactory() override = default;

  std::shared_ptr<const GenomeVariant> compoundSNP(const std::shared_ptr<const GenomeVariant>& SNPs,
                                                   const std::shared_ptr<const GenomeDatabase>& genome_db);

private:

  bool aggregateVariant(const std::shared_ptr<const Variant>& variant_ptr) const override;

};



}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_VARIANT_FACTORY_COMPOUND_H
