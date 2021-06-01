/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_GOA_ANNOTATION_PARSER
#define KGL_GOA_ANNOTATION_PARSER

#include "kol_ParserAnnotationInterface.h"
#include "kol_PolicyEvidence.h"



namespace kellerberrin::ontology {


/*! \class ParserAnnotationGoa
	\brief A class to parse a Uniprot Gene Ontolog Annotation (GOA) file.

	This class will read a GOA file an return an AnnotationData object pointer.
	  Defined at: http://www.ebi.ac.uk/GOA

	 Implements ParserAnnotationInterface

*/
class ParserAnnotationGoa : public ParserAnnotationInterface {

public:

  //! A default constructor method for creating the parser with a policy.
  /*!
    Creates the parser
  */
  ParserAnnotationGoa() = default;

  //! A parameterized constructor method for creating the parser with a policy.
  /*!
    Creates the parser
  */
  explicit ParserAnnotationGoa(const PolicyEvidence &policy) : _policy(policy) {}
  ~ParserAnnotationGoa() override = default;


  //! An interface method for parsing an annotation file.
  /*!
    This method takes a filename as in put and returns a pointer to an
      AnnotationData object. This method fulfills part of the interface contract.
  */
  [[nodiscard]] std::shared_ptr<AnnotationData> parseAnnotationFile(const std::string &filename) const override;

  //! A method for checking if a file exists and is formatted correctly.
  /*!
    This function checks that the file exists and its format can be recognized.
  */
  [[nodiscard]] bool isFileGood(const std::string &fileName) const override;

private:

  const PolicyEvidence _policy;

  const static constexpr char TAB_FIELD_DELIMITER = '\t';
  const static constexpr size_t EXPECTED_FIELD_COUNT = 17;
  const static constexpr char COMMENT_CHAR = '!';
  const static constexpr size_t FIELD_OFFSET_DB_ID = 0;
  const static constexpr size_t FIELD_OFFSET_UNIPROT_ID = 1;
  const static constexpr size_t FIELD_OFFSET_GENE_ID = 2;
  const static constexpr size_t FIELD_OFFSET_QUALIFIER = 3;
  const static constexpr char* NOT_QUALIFIER = "NOT";
  const static constexpr size_t FIELD_OFFSET_GO_ID = 4;
  const static constexpr size_t FIELD_OFFSET_DB_REFERENCE = 5;
  const static constexpr size_t FIELD_OFFSET_EVIDENCE_CODE = 6;
  const static constexpr size_t FIELD_OFFSET_WITH_FROM = 7;
  const static constexpr size_t FIELD_OFFSET_ONTOLOGY = 8;
  const static constexpr char* ONTOLOGY_BP_CODE = "P";
  const static constexpr char* ONTOLOGY_MF_CODE = "F";
  const static constexpr char* ONTOLOGY_CC_CODE = "C";
  const static constexpr size_t FIELD_OFFSET_NAME = 9;
  const static constexpr size_t FIELD_OFFSET_SYNONYM = 10;
  const static constexpr size_t FIELD_OFFSET_TYPE = 11;
  const static constexpr size_t FIELD_OFFSET_TAXON = 12;
  const static constexpr size_t FIELD_OFFSET_DATE = 13;
  const static constexpr size_t FIELD_OFFSET_ASSIGNED = 14;
  const static constexpr size_t FIELD_OFFSET_EXTENSION = 15;
  const static constexpr size_t FIELD_OFFSET_PRODUCT = 16;


};

} // namespace

#endif
