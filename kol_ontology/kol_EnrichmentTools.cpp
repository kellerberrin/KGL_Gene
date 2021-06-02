//
// Created by kellerberrin on 19/4/21.
//


#include "kol_OntologyTypes.h"
#include "kol_EnrichmentTools.h"

namespace kol = kellerberrin::ontology;

//! A method for determining which genes are annotated with the given term or a child of that term.
/*!
  This method calculates the set of the genes annotated with a given term or transatively with a child of that term.
*/
kol::OntologySetType<std::string> kol::EnrichmentTools::getDescendantGenes( const GoGraph &go,
                                                                            const TermAnnotation &data,
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
double kol::EnrichmentTools::enrichmentSignificance( const GoGraph &go,
                                                     const TermAnnotation &data,
                                                     const OntologySetType<std::string> &genes,
                                                     const std::string &term) {

  OntologySetType<std::string> termGenes = getDescendantGenes(go, data, term);
  OntologySetType<std::string> sharedGenes = SetUtilities::setIntersection(genes, termGenes);

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
