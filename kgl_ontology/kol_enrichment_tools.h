/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.

==============================================================================*/

#ifndef ENRICHMENT_TOOLS
#define ENRICHMENT_TOOLS

#include <stdlib.h>

#include "kel_distribution.h"
#include "kol_AnnotationData.h"
#include "kol_GoGraph.h"
#include "SetUtilities.h"

namespace kellerberrin::ontology {
  //! The EnrichmentTools namespace provides simple functions for calulating GO term enrichment.
  /*!
    This namespace defines free functions that allow enrichment p-values to be calculated.
    These funcitons can serve as the foundation for more sophisticated enrichment analysis.
  */
  class EnrichmentTools {

  public:

    // Just static functions.
    EnrichmentTools() = delete;

    //! A method for determining which genes are annotated with the given term or a child of that term.
    /*!
      This method calculates the set of the genes annotated with a given term or transatively with a child of that term.
    */
    [[nodiscard]] static OntologySetType<std::string> getDescendantGenes(const GoGraph &go,
                                                                         const AnnotationData &data,
                                                                         const std::string &term) {

      OntologySetType<std::string> descendants = go.getSelfDescendantTerms(term);

      OntologySetType<std::string> genes;
      for (auto const &current_term : descendants) {
        //std::cout << *si << " " << go->getTermName(*si) << std::endl;
        data.addGenesForGoTerm(current_term, genes);

      }

      return genes;

    }

    //! A method for calculating the result of a hypergeometic test.
    /*!
      This method calculates p-value of a hypergeometic test give 4 values.
      The number of successes_r successes_r in the total sample population_N where sample_size_n is the sample size without replacement.
      The population_N sample_size_n   successes_r
      The population_N size  population_N
      The test value           k

      Answers the question:
      "What is probability of seeing value of test_value_k or more successes
        in a sample of size sample_size_n, given that total sample size population_N contains total successes_r."
    */

    //! A method to calculate the enrichment of a specific term in a sample of genes.
    /*!
      This method performs a hypergeometic test of enrichment for a term given
      a set of genes that serves as the sample. The population is taken as all genes
      in the annotation database.
    */
    [[nodiscard]] [[maybe_unused]] static double enrichmentSignificance(const GoGraph &go,
                                                                        const AnnotationData &data,
                                                                        const OntologySetType<std::string> &genes,
                                                                        const std::string &term) {

      OntologySetType<std::string> termGenes = getDescendantGenes(go, data, term);
      OntologySetType<std::string> sharedGenes = SetUtilities::set_intersection(genes, termGenes);

      if (sharedGenes.empty()) {

        return 1.0;

      }

      size_t sampleSize_n = genes.size();
      size_t sampleWithTerm_k = sharedGenes.size();
      size_t populationWithTerm_K = termGenes.size();
      size_t populationSize_N = data.getNumGenes();

      HypergeometricDistribution hypergeometric(populationWithTerm_K, sampleSize_n, populationSize_N);
      return hypergeometric.upperSingleTailTest(sampleWithTerm_k);

    }

  };

} // namespace

#endif
