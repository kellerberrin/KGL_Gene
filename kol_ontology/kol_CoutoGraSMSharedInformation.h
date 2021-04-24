/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_GRASM_SHARED_INFORMATION
#define KGL_GRASM_SHARED_INFORMATION

#include "kol_SharedInformationInterface.h"
#include "kol_TermInformationContentMap.h"
#include "kol_GoGraph.h"
#include "kol_Accumulators.h"
#include "kol_SetUtilities.h"

#include <utility>
#include <algorithm>

#include <boost/graph/breadth_first_search.hpp>


namespace kellerberrin::ontology {

/*! \class CoutoGraSMSharedInformation
	\brief A class to calculate shared information across disjoint common ancestors using the exact algorithm as written in the paper.

	This class calculates shared information across disjoint common ancsetors.

    F. M. Couto, M. J. Silva, and P. M. Coutinho, "Measuring semantic similarity
	between Gene Ontology terms," Data & Knowledge Engineering, vol. 61, 
	pp. 137-152, Apr 2007.

	Couto proposed calculating this value as a subsitute for the IC of the MICA in calculating
	 Resnik, Lin, and Jiang-Conrath

*/
class CoutoGraSMSharedInformation : public SharedInformationInterface {

public:

  //! A constructor
  /*!
    Creates the CoutoGraSMGreaterOrEqual class
  */
  CoutoGraSMSharedInformation( const std::shared_ptr<const GoGraph> &graph_ptr,
                               const std::shared_ptr<const TermInformationContentMap> &ic_map_ptr)
      : graph_ptr_(graph_ptr), ic_map_ptr_(ic_map_ptr) {}

  ~CoutoGraSMSharedInformation() override = default;

  //! A method for determining the common disjunctive ancestors
  /*!
    This method returns the common disjunctive ancestors for two terms
  */
  [[nodiscard]] OntologySetType<std::string> getCommonDisjointAncestors(const std::string &termC1,
                                                                        const std::string &termC2) const;

  //! A method for determining if for a term c, a pair (a1,a2) is disjoint in c
  /*!
    This method returns
  */
  [[nodiscard]] bool isDisjoint(const std::string &termC, const std::string &termA1, const std::string &termA2) const;

  //! A method for calculating the number of paths for one term to another.
  /*!
    This method returns the number of paths between two terms
  */
  [[nodiscard]] size_t getNumPaths(const std::string &termA, const std::string &termB) const;


  //! An method for returning the shared information of two terms
  /*!
    This method returns the mean information content disjoint common ancestors
  */
  [[nodiscard]] double sharedInformation(const std::string &termA, const std::string &termB) const override;

  //! An interface method for returning the shared information of a single terms,or information content
  /*!
    This method privdes a mechanism for returing a term's information content.
  */
  [[nodiscard]] double sharedInformation(const std::string &term) const override;

  //! An interface method for returning the maximum information content for a term
  /*!
    This method provides the absolute max information content within a corpus for normalization purposes.
  */
  [[nodiscard]] double maxInformationContent(const std::string &term) const override;

  //! An interface method for determining if a term can be found
  /*!
    Determines if the term can be found in the current map.
  */
  [[nodiscard]] bool hasTerm(const std::string &term) const override {

    return ic_map_ptr_->hasTerm(term);

  }

  //! An interface method for determining if the two terms are of like ontologies.
  /*!
    Determine if two terms are of the same ontology.
  */
  [[nodiscard]] bool isSameOntology(const std::string &termA, const std::string &termB) const override {

    return graph_ptr_->getTermOntology(termA) == graph_ptr_->getTermOntology(termB);

  }


private:

  std::shared_ptr<const GoGraph> graph_ptr_;
  std::shared_ptr<const TermInformationContentMap> ic_map_ptr_;
  OntologyMapType<std::string, size_t> path_memory_;

  //! Count paths from B to A
  /*!
    Count paths between B and A
  */
  [[nodiscard]] size_t pathCount(const std::string &termA, const std::string &termB) const;

  //! Recursive helper method that performs the DFS topological sort for path counting
  /*!
    A path counting topological sort recursive method.
  */
  void visitHelper(const GoGraph::GoVertex &go_vertex,
                   const GoGraph::Graph &go_graph,
                   OntologySetType<std::string> &ancestors,
                   OntologySetType<std::string> &finished,
                   OntologyMapType<std::string, size_t> &pathMap) const;

  //! A private function to create a string key from a pair of terms
  /*!
    Creates a string key our of a pair to use in memorizing path counts
  */
  [[nodiscard]] std::string keyPair(const std::string &termA, const std::string &termB) const;

  //! A private function to test if the key as been seen already
  /*!
    A private function to test if the key as been seen already.
  */
  [[nodiscard]] bool hasSeenKey(const std::string &key) const;


};


} // namespace


#endif
