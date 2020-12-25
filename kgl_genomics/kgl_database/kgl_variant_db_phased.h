//
// Created by kellerberrin on 8/01/18.
//

#ifndef KGL_VARIANT_DB_PHASED_H
#define KGL_VARIANT_DB_PHASED_H




#include <map>
#include <memory>
#include <vector>
#include <sstream>
#include "kgl_attributes.h"
#include "kgl_variant.h"
#include "kgl_variant_mutation.h"
#include "kgl_variant_db_population.h"


namespace kellerberrin::genome {   //  organization level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// An Unphased Population. A General Population object for gnomad, clinvar etc.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class UnphasedBase : public DataObjectBase {

public:

  virtual ~UnphasedBase() override = default;

  DataTypeEnum dataType() const override { return DataTypeEnum::UnphasedPopulation; }

protected:

  explicit UnphasedBase(const PopulationId_t& population_id) : DataObjectBase(population_id) {}

};

// General unphased population.
using ContigVariant = ContigOffsetVariant;
using DiploidContig = ContigOffsetVariant;
using GenomeVariant = ContigOffsetVariant;
using UnphasedPopulation = PopulationVariant<UnphasedBase>;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A Haploid Population (P. Falciparum)
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class HaploidBase : public DataObjectBase {

public:

  virtual ~HaploidBase() override = default;

  DataTypeEnum dataType() const override { return DataTypeEnum::HaploidPopulation; }

protected:

  explicit HaploidBase(const PopulationId_t& population_id) : DataObjectBase(population_id) {}

};


// Haploid.
using HaploidPopulation = PopulationVariant<HaploidBase>;

// Used by mutation analytics.
// todo: Add templates to the mutation analytics classes
using PhasedPopulation = HaploidPopulation;
using PhasedGenome = GenomeVariantArray;
using HaploidGenome = GenomeVariantArray;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A Diploid Phased Population
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class DiploidBase : public DataObjectBase {

public:

  virtual ~DiploidBase() override = default;

  DataTypeEnum dataType() const override { return DataTypeEnum::DiploidPopulation; }

protected:

  explicit DiploidBase(const PopulationId_t& population_id) : DataObjectBase(population_id) {}

};


// Diploid.
using DiploidPopulation = PopulationVariant<DiploidBase>;


}   // end namespace



#endif //KGL_VARIANT_DB_PHASED_H
