/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KOL_GAF_ANNOTATION_PARSER
#define KOL_GAF_ANNOTATION_PARSER

#include "kol_PolicyEvidence.h"
#include "kol_TermAnnotation.h"

#include <memory>

namespace kellerberrin::ontology {


/*! \class ParserAnnotationGaf
	\brief A class to parse a Uniprot Gene Ontolog Annotation (GOA) file.

	This class will read a GOA file an return an TermAnnotation object pointer.
	  Defined at: http://www.ebi.ac.uk/GOA

	 Implements ParserAnnotationInterface

*/
class ParserAnnotationGaf  {

public:

  ParserAnnotationGaf() = delete;
  ~ParserAnnotationGaf() = delete;


  //! An interface method for parsing an annotation file.
  /*!
    This method takes a filename as in put and returns a pointer to an
      TermAnnotation object. This method fulfills part of the interface contract.
  */
  [[nodiscard]] static std::shared_ptr<const TermAnnotation> parseAnnotationFile( const PolicyEvidence &policy,
                                                                                  const std::string &filename,
                                                                                  AnnotationGeneName gene_id_type = AnnotationGeneName::UNIPROT_GENE_ID);

private:

  const static constexpr char TAB_FIELD_DELIMITER_ = '\t';
  const static constexpr size_t EXPECTED_FIELD_COUNT_ = 17;
  const static constexpr char COMMENT_CHAR_ = '!';
  const static constexpr size_t FIELD_OFFSET_QUALIFIER_ = 3;
  const static constexpr char* NOT_QUALIFIER_ = "NOT";

};

} // namespace

#endif
