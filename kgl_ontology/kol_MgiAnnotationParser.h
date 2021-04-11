/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef MGI_ANNOTATION_PARSER
#define MGI_ANNOTATION_PARSER

#include "kol_AnnotationParserInterface.h"

#include <iostream>
#include <boost/tokenizer.hpp>


namespace kellerberrin::ontology {

/*! \class MgiAnnotationParser
	\brief A class to parse an Mouse Genome Informatics go annotation file.

	This class will read an mgi annotation file and add those annoations to 
	  an AnnotationData class.
	  Available at: ftp://ftp.informatics.jax.org/pub/reports/index.html#go

	  MGI uses GAF format which is GOA.

	 Implements AnnotationParserInterface

*/
class MgiAnnotationParser : public GoaAnnotationParser {

public:
  //! A default constructor method for creating the parser
  /*!
    Creates the parser
  */
  MgiAnnotationParser() = default;

  //! A parameterized constructor method for creating the parser with a policy.
  /*!
    Creates the parser with a custom policy
  */
  explicit MgiAnnotationParser(const EvidencePolicyInterface &policy) : GoaAnnotationParser(policy) {}

  ~MgiAnnotationParser() override = default;

  [[nodiscard]] std::unique_ptr<AnnotationParserInterface> clone() const override { return std::make_unique<MgiAnnotationParser>(getPolicy()); }


};

}  // namespace

#endif
