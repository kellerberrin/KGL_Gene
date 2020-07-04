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



// General unphased population.
using ContigVariant = ContigOffsetVariant<UnphasedContigListOffset>;
using GenomeVariant = GenomeVariantArray<ContigVariant>;
using UnphasedPopulation = PopulationVariant<GenomeVariant>;

// Haploid.
using HaploidContig = ContigOffsetVariant<HaploidOffset>;
using HaploidGenome = GenomeVariantArray<HaploidContig>;
using HaploidPopulation = PopulationVariant<HaploidGenome>;

// Used by mutation analytics.
// todo: Add templates to the mutation analytics classes
using PhasedPopulation = HaploidPopulation;
using PhasedGenome = HaploidGenome;


// Diploid.
using DiploidContig = ContigOffsetVariant<DiploidOffset>;
using DiploidGenome = GenomeVariantArray<DiploidContig>;
using DiploidPopulation = PopulationVariant<DiploidGenome>;


}   // end namespace



#endif //KGL_VARIANT_DB_POPULATION_H
