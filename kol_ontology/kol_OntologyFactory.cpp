//
// Created by kellerberrin on 14/6/21.
//

#include "kol_OntologyFactory.h"
#include "kel_exec_env.h"

// Information
#include "kol_InformationExclusiveInherited.h"
#include "kol_InformationContentDAG.h"
#include "kol_InformationContent.h"
#include "kol_InformationCoutoGraSM.h"
#include "kol_InformationCoutoGraSMAdjusted.h"
#include "kol_InformationFrontier.h"
#include "kol_InformationAncestorMean.h"

// Similarity
#include "kol_SimilarityJiangConrath.h"
#include "kol_SimilarityResnik.h"
#include "kol_SimilarityRelevance.h"
#include "kol_SimilarityLin.h"
#include "kol_SimilarityPekarStaab.h"

// Set Similarity
#include "kol_SetSimilarityGentlemanSimUI.h"
#include "kol_SetSimilarityJaccard.h"
#include "kol_SetSimilarityMazanduSimDIC.h"
#include "kol_SetSimilarityMazanduSimUIC.h"
#include "kol_SetSimilarityPesquitaSimGIC.h"
#include "kol_SetSimilarityAverageBestMatch.h"
#include "kol_SetSimilarityBestMatchAverage.h"
#include "kol_SetSimilarityAllPairsMax.h"
#include "kol_SetSimilarityAllPairsAverage.h"



namespace kol = kellerberrin::ontology;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Implementation of the Information Factory.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::shared_ptr <const kol::InformationInterface>
    kol::OntologyFactory::createInformation( const std::shared_ptr<const TermAnnotation>& annotation_ptr,
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


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Implementation of the Similarity Factory.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



std::shared_ptr<const kol::SimilarityInterface>
    kol::OntologyFactory::createSimilarity( const std::shared_ptr<const TermAnnotation>& annotation_ptr,
                                            const std::shared_ptr<const GoGraph>& graph_ptr,
                                            SimilarityType similarity_type,
                                            InformationContentType information_type) {

  if (similarity_type == SimilarityType::PEKARSTAAB) {

    auto depth_info_ptr = std::make_shared<const InformationDepthMap>(*graph_ptr);
    return std::make_shared<const SimilarityPekarStaab>(graph_ptr, depth_info_ptr);

  }

  auto information_ptr = OntologyFactory::createInformation(annotation_ptr, graph_ptr, information_type);

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


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Implementation of the Set Similarity Factory.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::shared_ptr<const kol::SetSimilarityInterface> kol::OntologyFactory::createSetSimilarity( const std::shared_ptr<const TermAnnotation>& annotation_ptr,
                                                                                              const std::shared_ptr<const GoGraph>& graph_ptr,
                                                                                              SetSimilarityType set_similarity_type,
                                                                                              SimilarityType similarity_type,
                                                                                              InformationContentType information_type) {

  switch (set_similarity_type) {

    // These set similarity measures are stand-alone
    case SetSimilarityType::GENTLEMANSIMUI:
      return std::make_shared<const GentlemanSimUISetSimilarity>(graph_ptr);

    case SetSimilarityType::JACCARD:
      return std::make_shared<const JaccardSetSimilarity>();

      // These set similarity measures require an information source to be specified.
    case SetSimilarityType::MAZANDUSIMDIC:
      return std::make_shared<const SetSimilarityMazanduSimDIC>(graph_ptr, OntologyFactory::createInformation(annotation_ptr, graph_ptr, information_type));

    case SetSimilarityType::MAZANDUSIMUIC:
      return std::make_shared<const SetSimilarityMazanduSimUIC>(graph_ptr, OntologyFactory::createInformation(annotation_ptr, graph_ptr, information_type));

    case SetSimilarityType::PESQUITASIMGIC:
      return std::make_shared<const SetSimilarityPesquitaSimGIC>(graph_ptr, OntologyFactory::createInformation(annotation_ptr, graph_ptr, information_type));

      // These set similarity measures require a similarity measure to be specified.
    case SetSimilarityType::AVERAGEBESTMATCH:
      return std::make_shared<const SetSimilarityAverageBestMatch>(OntologyFactory::createSimilarity(annotation_ptr, graph_ptr, similarity_type, information_type));

    case SetSimilarityType::ALLPAIRSMAX:
      return std::make_shared<const SetSimilarityAllPairsMax>(OntologyFactory::createSimilarity(annotation_ptr, graph_ptr, similarity_type, information_type));

    case SetSimilarityType::ALLPAIRSAVERAGE:
      return std::make_shared<const SetSimilarityAllPairsAverage>(OntologyFactory::createSimilarity(annotation_ptr, graph_ptr, similarity_type, information_type));

    default:
    case SetSimilarityType::BESTMATCHAVERAGE:
      return std::make_shared<const SetSimilarityBestMatchAverage>(OntologyFactory::createSimilarity(annotation_ptr, graph_ptr, similarity_type, information_type));

  }

}




std::shared_ptr<const kol::SetSimilarityInterface> kol::OntologyFactory::createSetSimilarity( const std::shared_ptr<const SimilarityInterface>& similarity_ptr,
                                                                                              SetSimilarityType set_similarity_type) {

  switch (set_similarity_type) {

      // These set similarity measures require a similarity measure to be specified.
    case SetSimilarityType::AVERAGEBESTMATCH:
      return std::make_shared<const SetSimilarityAverageBestMatch>(similarity_ptr);

    case SetSimilarityType::ALLPAIRSMAX:
      return std::make_shared<const SetSimilarityAllPairsMax>(similarity_ptr);

    case SetSimilarityType::ALLPAIRSAVERAGE:
      return std::make_shared<const SetSimilarityAllPairsAverage>(similarity_ptr);

    case SetSimilarityType::BESTMATCHAVERAGE:
      return std::make_shared<const SetSimilarityBestMatchAverage>(similarity_ptr);

    default:
      ExecEnv::log().error("OntologyFactory::createSetSimilarity; Non SimilarityInterface type requested, BMA returned");
      return std::make_shared<const SetSimilarityBestMatchAverage>(similarity_ptr);

  }

}
