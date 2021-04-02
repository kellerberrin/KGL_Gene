/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef ENRICHMENT_TOOLS
#define ENRICHMENT_TOOLS

#include <stdlib.h>

#include <AnnotationData.h>
#include <GoGraph.h>
#include <SetUtilities.h>

#include <boost/math/distributions.hpp>
#include <boost/unordered_map.hpp>

//! The EnrichmentTools namespace provides simple functions for calulating GO term enrichment.
/*!
	This namespace defines free functions that allow enrichment p-values to be calculated.
	These funcitons can serve as the foundation for more sophisticated enrichment anlayis.
*/
class EnrichmentTools {

public:

  // Just static functions.
  EnrichmentTools() = delete;

	//! A method for determining which genes are annotated with the given term or a child of that term.
	/*!
		This method calculates the set of the genes annotated with a given term or transatively with a child of that term.
	*/
	[[nodiscard]] static OntologySetType<std::string> getDescendantGenes( const GoGraph &go,
                                                                        const AnnotationData& data,
                                                                        const std::string &term) {

		OntologySetType<std::string> descendants = go.getDescendantTerms(term);
		descendants.insert(term);

    OntologySetType<std::string> genes;

		for(auto const& current_term : descendants) {
			//std::cout << *si << " " << go->getTermName(*si) << std::endl;
			data.addGenesForGoTerm(current_term, genes);

		}
		//std::cout << genes.size() << std::endl;

		return genes;

	}

	//! A method for calculating the result of a hypergeometic test.
	/*!
		This method calculates p-value of a hypergeometice test give 4 values.
		The sample size,         n
		The population success   K
		The the population size  N
		The test value           k

		Answers the question:
		"What is probability of seeing value of k or more successes
		  in a sample of size n, given that the population of size N 
		  contains K total successes."
	*/
  [[nodiscard]] static double oneSidedRawPvalue_hyper(size_t sample, size_t success, size_t population, size_t test_value) {

		double sum = 0.0;
		boost::math::hypergeometric dist(sample,success,population);
		for(size_t i = test_value; i <= sample && i <= success; ++i){

			double prob = boost::math::pdf(dist,i);
			sum += prob;

		}

		return sum;

	}


	//! A method to calculate the enrichment of a specific term in a sample of genes.
	/*!
		This method performs a hypergeometic test of enrichment for a term given
		a set of genes that serves as the sample. The population is taken as all genes
		in the annotation database.
	*/
	[[nodiscard]] [[maybe_unused]] static double enrichmentSignificance( const GoGraph& go,
                                                                       const AnnotationData& data,
                                                                       const OntologySetType<std::string>& genes,
											                                                 const std::string &term)
	{

		OntologySetType<std::string> termGenes = getDescendantGenes(go, data, term);
    OntologySetType<std::string> sharedGenes = SetUtilities::set_intersection(genes,termGenes);

		if(sharedGenes.empty()){

			return 1.0;

		}
		
		size_t sampleSize = genes.size();
		size_t sampleWithTerm = sharedGenes.size();
		size_t populationWithTerm = termGenes.size();
		size_t populationSize = data.getNumGenes();

		return oneSidedRawPvalue_hyper(sampleSize,populationWithTerm,populationSize,sampleWithTerm);

	}

};

#endif
