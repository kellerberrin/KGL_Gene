/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef PEKAR_STAAB_SIMILARITY
#define PEKAR_STAAB_SIMILARITY

#include <kol_TermSimilarityInterface.h>
#include <kol_TermDepthMap.h>
#include <kol_GoGraph.h>

namespace kellerberrin::ontology {

/*! \class PekarStaabSimilarity
	\brief A class to calculate PekarStaab similarity between 2 terms

	This class calculates Pekar Staab similarity.
	
	V. Pekar and S. Staab, "Taxonomy learning: factoring the structure 
	 of a taxonomy into a semantic classification decision," in 
	 Proc. of 19th International Conference on Computational Linguistics. 
	 Morristown NJ USA: Association for Computational Linguistics, pp. 1-7, 2002.

	H. Yu, L. Gao, K. Tu, and Z. Guo, "Broadly predicting specific gene 
	 functions with expression similarity and taxonomy similarity,"
	 Gene, vol. 352, pp. 75-81, Jun 6 2005.
	  
    lowest common ancestor (LCA)
	GraphDist(LCA,root)/(GraphDist(a,LCA)+GraphDist(b,LCA)+GraphDist(LCA,root))

*/
class PekarStaabSimilarity : public TermSimilarityInterface {

public:

  //! A constructor
  /*!
    Creates the default(empty) StandardRelationshipPolicy
  */
  PekarStaabSimilarity( const std::shared_ptr<const GoGraph> &graph_ptr,
                        const std::shared_ptr<const TermDepthMap> &depth_map_ptr)
      : graph_ptr_(graph_ptr), depth_map_ptr_(depth_map_ptr) {}

  ~PekarStaabSimilarity() override = default;


  //! A method for calculating term-to-term similarity for GO terms using Pekar Staab similarity
  /*!
    This method returns the PekarStaab similarity.
  */
  [[nodiscard]] double calculateTermSimilarity(const std::string &goTermA, const std::string &goTermB) const override;

  //! A method for calculating term-to-term similarity for GO terms using Normalized Pekar Staab similarity
  /*!
    This method returns the PekarStaab similarity scaled between 0 and 1 [0,1] inclusive
  */
  [[nodiscard]] double calculateNormalizedTermSimilarity(const std::string &goTermA, const std::string &goTermB) const override {
    //Pekar and Staab's method is already normalized
    return calculateTermSimilarity(goTermA, goTermB);

  }


private:

  std::shared_ptr<const GoGraph> graph_ptr_;
  std::shared_ptr<const TermDepthMap> depth_map_ptr_;

};


} // namespace


#endif