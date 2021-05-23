/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_GAF_ANNOTATION_PARSER
#define KGL_GAF_ANNOTATION_PARSER

#include "kol_ParserAnnotationGoa.h"

namespace kellerberrin::ontology {


/*! \class ParserAnnotationGaf
\brief A class to parse a GO Annotation File (GAF, Format 2.0).

This class will read a GAF file and return an AnnotationData object pointer.
Defined at: http://geneontology.org/page/go-annotation-file-format-20

Implements ParserAnnotationInterface.

For now, the important aspects of the GAF file and GOA file are the same.
The ParserAnnotationGaf inherits all functionality from ParserAnnotationGoa. This may change in the future.
*/
class ParserAnnotationGaf : public ParserAnnotationGoa {

public:
  //! A default constructor method for creating the parser
  /*!
    Creates the parser
  */
  ParserAnnotationGaf() : ParserAnnotationGoa() {}

  //! A parameterized constructor method for creating the parser with a policy.
  /*!
    Creates the parser with a custom policy
  */
  explicit ParserAnnotationGaf(const PolicyAllowedEvidence &policy) : ParserAnnotationGoa(policy) {}
  ~ParserAnnotationGaf() override = default;


};

}  // namespace

#endif
