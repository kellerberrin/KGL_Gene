//
// Created by kellerberrin on 10/6/21.
//

#ifndef KOL_INFORMATIONFACTORY_H
#define KOL_INFORMATIONFACTORY_H


#include "kol_InformationExclusiveInherited.h"
#include "kol_InformationContentDAG.h"
#include "kol_InformationContent.h"
#include "kol_InformationCoutoGraSM.h"
#include "kol_InformationCoutoGraSMAdjusted.h"
#include "kol_InformationFrontier.h"
#include "kol_InformationAncestorMean.h"



namespace kellerberrin::ontology {


enum class InformationContentType {
  ANCESTOR, CONTENT, CONTENT_DAG, COUTOGRASM, COUTOGRASMADJUSTED, EXCLUSIVEINHERITED, FRONTIER
};


class InformationContentFactory {

public:

  /*!
    This object cannot be created.
  */
  InformationContentFactory() = delete;
  ~InformationContentFactory() = delete;


  [[nodiscard]] static std::shared_ptr <const InformationInterface> createInformation( const std::shared_ptr<const TermAnnotation>& annotation_ptr,
                                                                                       const std::shared_ptr<const GoGraph>& graph_ptr,
                                                                                       InformationContentType information_type) {

    switch (information_type) {

      case InformationContentType::CONTENT_DAG:
        return std::make_shared<const InformationContentDAG>(graph_ptr, annotation_ptr);

      case InformationContentType::ANCESTOR: {

        auto info_content_ptr = std::make_shared<const InformationContent>(graph_ptr, annotation_ptr);
        return std::make_shared<const InformationAncestorMean>(graph_ptr, info_content_ptr);

      }

      case InformationContentType::COUTOGRASM: {

        auto info_content_ptr = std::make_shared<const InformationContent>(graph_ptr, annotation_ptr);
        return std::make_shared<const InformationCoutoGraSM>(graph_ptr, info_content_ptr);

      }

      case InformationContentType::COUTOGRASMADJUSTED: {

        auto info_content_ptr = std::make_shared<const InformationContent>(graph_ptr, annotation_ptr);
        return std::make_shared<const InformationCoutoGraSMAdjusted>(graph_ptr, info_content_ptr);

      }

      case InformationContentType::EXCLUSIVEINHERITED: {

        auto info_content_ptr = std::make_shared<const InformationContent>(graph_ptr, annotation_ptr);
        return std::make_shared<const InformationExclusiveInherited>(graph_ptr, info_content_ptr);

      }

      case InformationContentType::FRONTIER: {

        auto info_content_ptr = std::make_shared<const InformationContent>(graph_ptr, annotation_ptr);
        return std::make_shared<const InformationFrontier>(graph_ptr, info_content_ptr);

      }

      default:
      case InformationContentType::CONTENT:
        return std::make_shared<const InformationContent>(graph_ptr, annotation_ptr);

    }

  }


};


} // namespace


#endif //KOL_INFORMATIONFACTORY_H
