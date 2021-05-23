/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_ANNOTATION_PARSER_FACTORY
#define KGL_ANNOTATION_PARSER_FACTORY

#include "kol_ParserAnnotationInterface.h"
#include "kol_ParserAnnotationEntrezGene2Go.h"
#include "kol_ParserAnnotationGoa.h"
#include "kol_ParserAnnotationGaf.h"
#include "kol_ParserAnnotationMgi.h"
#include <vector>
#include <string>


namespace kellerberrin::ontology {

/*! \class ParserAnnotationFactory
	\brief A class to return an instance of ParserAnnotationInterface at runtime based on an argument.

	This class holds a set of parser classes. When queried using the getParser method, it
	  returns an instance of ParserAnnotationInterface based on a string key. This allows
	  parsers to be easily added to a larger system and switched at runtime.

*/
enum class AnnotationParserType { ENTREZ_ANNO_PARSER, GOA_ANNO_PARSER, GAF_ANNO_PARSER, MGI_ANNO_PARSER };

class ParserAnnotationFactory {


public:

  // Cannot create object.
  ParserAnnotationFactory() = delete;
  ~ParserAnnotationFactory() = delete;


  [[nodiscard]] static std::unique_ptr<ParserAnnotationInterface> createAnnotationParser(AnnotationParserType parser_type,
                                                                                         const PolicyAllowedEvidence &policy = PolicyAllowedEvidence()) {

    switch (parser_type) {

      case AnnotationParserType::ENTREZ_ANNO_PARSER:
        return std::make_unique<ParserAnnotationEntrezGene2Go>(policy);

      case AnnotationParserType::GOA_ANNO_PARSER:
        return std::make_unique<ParserAnnotationGoa>(policy);

      case AnnotationParserType::MGI_ANNO_PARSER:
        return std::make_unique<ParserAnnotationMgi>(policy);

      default:
      case AnnotationParserType::GAF_ANNO_PARSER:
        return std::make_unique<ParserAnnotationGaf>(policy);

    }

  }

};

} // namespace

#endif