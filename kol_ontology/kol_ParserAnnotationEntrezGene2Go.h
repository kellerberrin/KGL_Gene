/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_ENTREZ_GENE2GO_ANNOTATION_PARSER
#define KGL_ENTREZ_GENE2GO_ANNOTATION_PARSER

#include "kol_ParserAnnotationInterface.h"
#include "kol_PolicyEvidence.h"


namespace kellerberrin::ontology {

/*! \class ParserAnnotationEntrezGene2Go
	\brief A class to parse an Entrez gene2go annotation file.

	This class will read an Entrez gene2go file and add those annoations to 
	  an AnnotationData class.
	  Available at: ftp://ftp.ncbi.nih.gov/gene/DATA/

	 Implements ParserAnnotationInterface

*/
class ParserAnnotationEntrezGene2Go : public ParserAnnotationInterface {

public:

  //! A default constructor method for creating the parser with the default policy.
  /*!
    Creates the parser with the default evidence policy, everything is allowed.
  */
  ParserAnnotationEntrezGene2Go() = default;

  //! A parameterized constructor method for creating the parser with a policy.
  /*!
    Creates the parser
  */
  ParserAnnotationEntrezGene2Go(const PolicyEvidence &policy) : policy_(policy) {}
  ~ParserAnnotationEntrezGene2Go() override = default;


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

  const PolicyEvidence policy_;

};

} // namespace

#endif
