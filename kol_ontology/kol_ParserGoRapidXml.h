/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_RAPIDXML_GO_PARSER
#define KGL_RAPIDXML_GO_PARSER

#include "kol_ParserGoInterface.h"

#include <vector>
#include <string>

#include <xml/rapidxml_utils.h>
#include <xml/rapidxml.h>

#include <boost/tokenizer.hpp>

namespace kellerberrin::ontology {

/*! \class GoParserRapidXml
	\brief This class parses a go XML file using RapidXML library.

	This class parses a Gene Ontology XML file using RapidXML library.

	Implements ParserGoInterface

*/
class GoParserRapidXml : public ParserGoInterface {

public:

  GoParserRapidXml() = default;

  ~GoParserRapidXml() override = default;



  //! Method to parse the go file, should be an XML file.
  /*!
    This method will read a Gene Ontology XML file and add all relationships
     to the graph.

  */

  [[nodiscard]] std::unique_ptr<GoGraph> parseGoFile(const std::string &filename) const override;

  //! A method to test if a file fits the accepted format
  /*!
  Returns true if the file matches accepted format, false otherwise
  */
  [[nodiscard]] bool isFileGood(const std::string &filename) const override;

private:


};


} // namespace


#endif
