//
// Created by kellerberrin on 10/6/21.
//

#ifndef KOL_SIMILARITYFACTORY_H
#define KOL_SIMILARITYFACTORY_H

#include "kol_InformationFactory.h"
#include "kol_SimilarityJiangConrath.h"
#include "kol_SimilarityResnik.h"
#include "kol_SimilarityRelevance.h"
#include "kol_SimilarityLin.h"
#include "kol_SimilarityPekarStaab.h"



namespace kellerberrin::ontology {


enum class SimilarityType { RESNIK, LIN, RELEVANCE, JIANG, PEKARSTAAB };


class SimilarityFactory {

public:

  /*!
    This object cannot be created.
  */
  SimilarityFactory() = delete;
  ~SimilarityFactory() = delete;


  [[nodiscard]] static std::shared_ptr <const SimilarityInterface> createSimilarity( const std::shared_ptr<const TermAnnotation>& annotation_ptr,
                                                                                     const std::shared_ptr<const GoGraph>& graph_ptr,
                                                                                     SimilarityType similarity_type
                                                                                     InformationContentType information_type == InformationContentType::CONTENT) {

    if (similarity_type == SimilarityType::PEKARSTAAB) {

      return createPekarStaab(graph_ptr);

    }

    auto information_ptr = InformationContentFactory::createInformation(graph_ptr, annotation_ptr, information_type);

    switch (similarity_type) {

      case SimilarityType::RESNIK:
        return std::make_shared<const SimilarityResnik>(information_ptr);

      case SimilarityType::LIN:
        return std::make_shared<const SimilarityLin>(information_ptr);

      case SimilarityType::RELEVANCE:
        return std::make_shared<const SimilarityRelevance>(information_ptr);

      default:
      case SimilarityType::JIANG:
        return std::make_shared<const SimilarityJiangConrath>(information_ptr);

    }

  }

private:

  [[nodiscard]] static std::shared_ptr <const SimilarityInterface> createPekarStaab( const std::shared_ptr<const GoGraph>& graph_ptr) {

    auto depth_info_ptr = std::make_shared<const InformationDepthMap>(*graph_ptr);
    return std::make_shared<const SimilarityPekarStaab>(graph_ptr, depth_info_ptr);

  }

};


} // namespace



#endif //KOL_SIMILARITYFACTORY_H
