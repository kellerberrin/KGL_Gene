/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_ENTREZ_GENE2GO_ANNOTATION_PARSER
#define KGL_ENTREZ_GENE2GO_ANNOTATION_PARSER

#include "kol_AnnotationParserInterface.h"
#include "kol_DisallowedSetEvidencePolicy.h"


namespace kellerberrin::ontology {

/*! \class EntrezGene2GoAnnotationParser
	\brief A class to parse an Entrez gene2go annotation file.

	This class will read an Entrez gene2go file and add those annoations to 
	  an AnnotationData class.
	  Available at: ftp://ftp.ncbi.nih.gov/gene/DATA/

	 Implements AnnotationParserInterface

*/
class EntrezGene2GoAnnotationParser : public AnnotationParserInterface {

public:

  //! A default constructor method for creating the parser with the default policy.
  /*!
    Creates the parser with the default evidence policy, everything is allowed.
  */
  EntrezGene2GoAnnotationParser() : policy_ptr_(std::make_unique<const DisallowedSetEvidencePolicy>()) {}

  //! A parameterized constructor method for creating the parser with a policy.
  /*!
    Creates the parser
  */
  EntrezGene2GoAnnotationParser(const EvidencePolicyInterface &policy) : policy_ptr_(policy.clone()) {}

  ~EntrezGene2GoAnnotationParser() override = default;

  //! An interface method for creating a new instance of the parser.
  /*!
    This method returns a new instance of the class. This method partially
      fulfills the interface contract.
  */
  [[nodiscard]] std::unique_ptr<AnnotationParserInterface> clone() const override { return std::make_unique<EntrezGene2GoAnnotationParser>(*policy_ptr_); }


  //! An interface method for parsing an annotation file.
  /*!
    This method takes a filename as in put and returns a pointer to an
      AnnotationData object. This method fulfills part of the interface contract.
  */
  [[nodiscard]] std::unique_ptr<AnnotationData> parseAnnotationFile(const std::string &filename) const override;

  //! A method for checking if a file exists and is formatted correctly.
  /*!
    This function checks that the file exists and its format can be recognized.
  */
  [[nodiscard]] bool isFileGood(const std::string &fileName) const override;

private:

  std::unique_ptr<const EvidencePolicyInterface> policy_ptr_;

};

} // namespace

#endif
