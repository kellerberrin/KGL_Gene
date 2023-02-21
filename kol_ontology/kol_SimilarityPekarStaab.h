/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef KOL_PEKAR_STAAB_SIMILARITY
#define KOL_PEKAR_STAAB_SIMILARITY

#include <kol_SimilarityInterface.h>
#include <kol_InformationDepthMap.h>
#include <kol_GoGraph.h>

namespace kellerberrin::ontology {

/*! \class SimilarityPekarStaab
	\brief A class to calculate PekarStaab similarity between 2 terms

	This class calculates Pekar Staab similarity.
	
	V. Pekar and S. Staab, "Taxonomy learning: factoring the structure 
	 of a taxonomy into a semantic classification decision," in 
	 ThreadFunc. of 19th International Conference on Computational Linguistics.
	 Morristown NJ USA: Association for Computational Linguistics, pp. 1-7, 2002.

	H. Yu, L. Gao, K. Tu, and Z. Guo, "Broadly predicting specific gene 
	 functions with expression similarity and taxonomy similarity,"
	 Gene, vol. 352, pp. 75-81, Jun 6 2005.
	  
    lowest common ancestor (LCA)
	GraphDist(LCA,root)/(GraphDist(a,LCA)+GraphDist(b,LCA)+GraphDist(LCA,root))

*/
class SimilarityPekarStaab : public SimilarityInterface {

public:

  //! A constructor
  /*!
    Creates the default(empty) PolicyStandardRelationship
  */
  SimilarityPekarStaab(const std::shared_ptr<const GoGraph> &graph_ptr,
                       const std::shared_ptr<const InformationDepthMap> &depth_map_ptr)
      : graph_ptr_(graph_ptr), depth_map_ptr_(depth_map_ptr) {}

  ~SimilarityPekarStaab() override = default;


  //! A method for calculating term-to-term similarity for GO terms using Pekar Staab similarity
  /*!
    This method returns the PekarStaab similarity.
  */
  [[nodiscard]] double calculateTermSimilarity(const std::string &goTermA, const std::string &goTermB) const override;


private:

  std::shared_ptr<const GoGraph> graph_ptr_;
  std::shared_ptr<const InformationDepthMap> depth_map_ptr_;

};


} // namespace


#endif