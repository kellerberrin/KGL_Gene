/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KOL_TERM_SIMILARITY_INTERFACE
#define KOL_TERM_SIMILARITY_INTERFACE

#include <string>
#include "kol_GoEnums.h"


namespace kellerberrin::ontology {

//! An interface class for comparing semantic similarity of GO terms
/*! \class SimilarityInterface

	This class defines the interface for comparing term-to-term GO similarity.
*/
class SimilarityInterface {

public:

  SimilarityInterface() = default;

  virtual ~SimilarityInterface() = default;

  //! A pure virtual method for calculating term-to-term similarity for GO terms
  /*!
    This pure virtual method requires similarity measure to implement the basic interface
      that returns a similarity value for two go terms.
  */
  [[nodiscard]] virtual double calculateTermSimilarity(const std::string &go_termA, const std::string &go_termB) const = 0;


};

} // namespace

#endif
