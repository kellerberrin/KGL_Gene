/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_ANNOTATION_PARSER_FACTORY
#define KGL_ANNOTATION_PARSER_FACTORY

#include "kol_AnnotationParserInterface.h"
#include "DisallowedSetEvidencePolicy.h"
#include "kol_EntrezGene2GoAnnotationParser.h"
#include "kol_GoaAnnotationParser.h"
#include "kol_GafAnnotationParser.h"
#include "kol_MgiAnnotationParser.h"
#include <vector>
#include <string>


namespace kellerberrin::ontology {

/*! \class AnnotationParserFactory
	\brief A class to return an instance of AnnotationParserInterface at runtime based on an argument.

	This class holds a set of parser classes. When queried using the getParser method, it
	  returns an instance of AnnotationParserInterface based on a string key. This allows
	  parsers to be easily added to a larger system and switched at runtime.

*/
enum class AnnotationParserType {
  ENTREZ_ANNO_PARSER, GOA_ANNO_PARSER, GAF_ANNO_PARSER, MGI_ANNO_PARSER
};

class AnnotationParserFactory {


public:

  // Cannot create object.
  AnnotationParserFactory() = delete;

  ~AnnotationParserFactory() = delete;


  [[nodiscard]] static std::unique_ptr<AnnotationParserInterface> createAnnotationParser(AnnotationParserType parser_type,
                                                                                         const EvidencePolicyInterface &policy = DisallowedSetEvidencePolicy()) {

    switch (parser_type) {

      case AnnotationParserType::ENTREZ_ANNO_PARSER:
        return std::make_unique<EntrezGene2GoAnnotationParser>(policy);

      case AnnotationParserType::GOA_ANNO_PARSER:
        return std::make_unique<GoaAnnotationParser>(policy);

      case AnnotationParserType::MGI_ANNO_PARSER:
        return std::make_unique<MgiAnnotationParser>(policy);

      default:
      case AnnotationParserType::GAF_ANNO_PARSER:
        return std::make_unique<GafAnnotationParser>(policy);

    }

  }

};

} // namespace

#endif