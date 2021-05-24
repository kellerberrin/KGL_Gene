/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef MGI_ANNOTATION_PARSER
#define MGI_ANNOTATION_PARSER

#include "kol_ParserAnnotationInterface.h"

#include <iostream>
#include <boost/tokenizer.hpp>


namespace kellerberrin::ontology {

/*! \class ParserAnnotationMgi
	\brief A class to parse an Mouse Genome Informatics go annotation file.

	This class will read an mgi annotation file and add those annotations to
	  an AnnotationData class.
	  Available at: ftp://ftp.informatics.jax.org/pub/reports/index.html#go

	  MGI uses GAF format which is GOA.

	 Implements ParserAnnotationInterface

*/
class ParserAnnotationMgi : public ParserAnnotationGoa {

public:
  //! A default constructor method for creating the parser
  /*!
    Creates the parser
  */
  ParserAnnotationMgi() = default;

  //! A parameterized constructor method for creating the parser with a policy.
  /*!
    Creates the parser with a custom policy
  */
  explicit ParserAnnotationMgi(const PolicyEvidence &policy) : ParserAnnotationGoa(policy) {}
  ~ParserAnnotationMgi() override = default;

};

}  // namespace

#endif
