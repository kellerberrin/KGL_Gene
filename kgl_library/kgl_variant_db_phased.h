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

  DataTypeEnum dataType() override { return DataTypeEnum::UnphasedPopulation; }

protected:

  explicit UnphasedBase(const PopulationId_t& population_id) : DataObjectBase(population_id) {}

};

// General unphased population.
using ContigVariant = ContigOffsetVariant<UnphasedContigListOffset>;
using GenomeVariant = GenomeVariantArray<ContigVariant>;
using UnphasedPopulation = PopulationVariant<GenomeVariant, UnphasedBase>;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A Haploid Population (P. Falciparum)
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class HaploidBase : public DataObjectBase {

public:

  virtual ~HaploidBase() override = default;

  DataTypeEnum dataType() override { return DataTypeEnum::HaploidPopulation; }

protected:

  explicit HaploidBase(const PopulationId_t& population_id) : DataObjectBase(population_id) {}

};


// Haploid.
using HaploidContig = ContigOffsetVariant<HaploidOffset>;
using HaploidGenome = GenomeVariantArray<HaploidContig>;
using HaploidPopulation = PopulationVariant<HaploidGenome, HaploidBase>;

// Used by mutation analytics.
// todo: Add templates to the mutation analytics classes
using PhasedPopulation = HaploidPopulation;
using PhasedGenome = HaploidGenome;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A Diploid Phased Population
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class DiploidBase : public DataObjectBase {

public:

  virtual ~DiploidBase() override = default;

  DataTypeEnum dataType() override { return DataTypeEnum::DiploidPopulation; }

protected:

  explicit DiploidBase(const PopulationId_t& population_id) : DataObjectBase(population_id) {}

};


// Diploid.
using DiploidContig = ContigOffsetVariant<DiploidOffset>;
using DiploidGenome = GenomeVariantArray<DiploidContig>;
using DiploidPopulation = PopulationVariant<DiploidGenome, DiploidBase>;


}   // end namespace



#endif //KGL_VARIANT_DB_PHASED_H
