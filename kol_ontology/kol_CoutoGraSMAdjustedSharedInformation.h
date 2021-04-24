/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_GRASM_SHARED_INFORMATION_ADJUSTED
#define KGL_GRASM_SHARED_INFORMATION_ADJUSTED

#include "kol_SharedInformationInterface.h"
#include "kol_TermInformationContentMap.h"
#include "kol_GoGraph.h"
#include "kol_SetUtilities.h"
#include "kol_Accumulators.h"

#include <utility>
#include <algorithm>

#include <boost/graph/breadth_first_search.hpp>

namespace kellerberrin::ontology {


/*! \class CoutoGraSMAdjustedSharedInformation
	\brief A class to calculate shared infromation accross disjoint common ancetors using an adjusted algorithm.

	This class calculates shared infromation accross disjoint common ancetors. This is a modification of the
	 original algorithm provided by Couto. The adjustment changes the contrain to path lengths to strictly greater than.
	 See line 150.

    F. M. Couto, M. J. Silva, and P. M. Coutinho, "Measuring semantic similarity
	between Gene Ontology terms," Data & Knowledge Engineering, vol. 61, 
	pp. 137-152, Apr 2007.

	Couto proposing calculating this value a subsituite for the IC of the MICA in calculating
	 Resnik, Lin, and Jiang-Conrath

*/
class CoutoGraSMAdjustedSharedInformation : public SharedInformationInterface {

public:

  //! Constructor
  /*!
    Creates the CoutoGraSMAdjustedSharedInformation class
  */
  CoutoGraSMAdjustedSharedInformation( const std::shared_ptr<const GoGraph> &graph_ptr,
                                       const std::shared_ptr<const TermInformationContentMap> &ic_map_ptr_)
      : graph_ptr_(graph_ptr), ic_map_ptr_(ic_map_ptr_) {}

  ~CoutoGraSMAdjustedSharedInformation() override = default;

  //! Calculate disjunctive ancestors.
  /*!
    A method for determining common disjunctive ancestors for two terms
  */
  [[nodiscard]] OntologySetType<std::string> getCommonDisjointAncestors( const std::string &termC1,
                                                                         const std::string &termC2) const;

  //! Determine if a terms are disjoint in a concept.
  /*!
    A method for determining if, for a term c, a pair (a1,a2) is disjoint in c
  */
  [[nodiscard]] bool isDisjoint(const std::string &termC, const std::string &termA1, const std::string &termA2) const;

  //! Calculate the number of paths between two concept terms.
  /*!
    A method for calculating the number of paths from one term to another.
  */
  [[nodiscard]] size_t getNumPaths(const std::string &termA, const std::string &termB) const;

  //! Shared infromation between two conecepts.
  /*!
    A method for calculating the shared infromation between two concepts.
  */
  [[nodiscard]] double sharedInformation(const std::string &termA, const std::string &termB) const override;

  //! Term information content.
  /*!
    An interface method to conventiently get information content of a single term
  */
  [[nodiscard]] double sharedInformation(const std::string &term) const override;

  //! Maximum Ontology IC for normalization.
  /*!
    An interface method for returning the maximum information content for a term within a corpus for normalization purposes.
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
  [[nodiscard]] std::size_t pathCount(const std::string &termA, const std::string &termB) const;

  //! Recursive helper method that performs the DFS topological sort for path counting
  /*!
    A path counting topological sort recursive method.
  */
  void visitHelper(const GoGraph::GoVertex &v,
                   const GoGraph::Graph &graph,
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
