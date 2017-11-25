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

  virtual std::shared_ptr<GenomeVariant> disaggregate(const std::shared_ptr<const GenomeVariant>& genome_variant,
                                                      const std::shared_ptr<const GenomeDatabase>& genome_db_ptr) const = 0;

  virtual std::shared_ptr<GenomeVariant> create(const std::shared_ptr<const GenomeVariant>& genome_variants,
                                                const std::shared_ptr<const GenomeDatabase>& genome_db_ptr) const = 0;


};





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The top level insert/delete compound variant factory.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VariantInsertDeleteFactory : public VariantCompoundFactory {

public:

  explicit VariantInsertDeleteFactory() = default;
  ~VariantInsertDeleteFactory() override = default;

  std::shared_ptr<GenomeVariant> create(const std::shared_ptr<const GenomeVariant>& genome_variants,
                                        const std::shared_ptr<const GenomeDatabase>& genome_db_ptr) const override;

private:


  bool aggregateVariants(const std::shared_ptr<const GenomeVariant>& variant_ptr,
                         std::vector<std::shared_ptr<const CompoundVariantMap>>& aggregated_variants_vec) const;

  // Determine which variants to aggregate.
  virtual bool selectVariant(const std::shared_ptr<const Variant>& variant_ptr) const = 0;

  virtual std::shared_ptr<const Variant> createCompoundVariant(const CompoundVariantMap& variant_map) const = 0;

};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The compound delete variant factory.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VariantDeleteFactory : public VariantInsertDeleteFactory {

public:

  explicit VariantDeleteFactory() = default;
  ~VariantDeleteFactory() override = default;


private:

  bool selectVariant(const std::shared_ptr<const Variant>& variant_ptr) const override;

  std::shared_ptr<const Variant> createCompoundVariant(const CompoundVariantMap& variant_map) const override;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The compound insert variant factory.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VariantInsertFactory : public VariantInsertDeleteFactory {

public:

  explicit VariantInsertFactory() = default;
  ~VariantInsertFactory() override = default;

private:

  bool selectVariant(const std::shared_ptr<const Variant>& variant_ptr) const override;

  std::shared_ptr<const Variant> createCompoundVariant(const CompoundVariantMap& variant_map) const override;


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The compound SNP variant factory.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VariantCompoundSNPFactory : public VariantCompoundFactory {

public:

  explicit VariantCompoundSNPFactory() = default;
  ~VariantCompoundSNPFactory() override = default;

  std::shared_ptr<GenomeVariant> create(const std::shared_ptr<const GenomeVariant>& genome_variants,
                                        const std::shared_ptr<const GenomeDatabase>& genome_db_ptr) const override;

private:


  bool aggregateVariants(const std::shared_ptr<const GenomeVariant>& variant_ptr,
                         std::vector<std::shared_ptr<const SNPCompoundVariantMap>>& aggregated_variants_vec) const;

  std::shared_ptr<const Variant> createCompoundVariant(const SNPCompoundVariantMap& variant_map) const;

};



}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_VARIANT_FACTORY_COMPOUND_H
