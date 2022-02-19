//
// Created by kellerberrin on 14/6/21.
//

#ifndef KOL_SIMILARITYFACTORYIMPL_H
#define KOL_SIMILARITYFACTORYIMPL_H


#include "kol_GoGraph.h"
#include "kol_TermAnnotation.h"
#include "kol_InformationInterface.h"
#include "kol_SimilarityInterface.h"
#include "kol_SetSimilarityInterface.h"



namespace kellerberrin::ontology {


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Information.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum class InformationContentType {
  ANCESTOR, CONTENT, CONTENT_DAG, COUTOGRASM, COUTOGRASMADJUSTED, EXCLUSIVEINHERITED, FRONTIER
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Similarity.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum class SimilarityType { RESNIK, LIN, RELEVANCE, JIANG, PEKARSTAAB };

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SetSimilarity.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


enum class SetSimilarityType { GENTLEMANSIMUI, JACCARD,   // Stand-alone.
  MAZANDUSIMDIC, MAZANDUSIMUIC, PESQUITASIMGIC,   // Information source.
  BESTMATCHAVERAGE, AVERAGEBESTMATCH, ALLPAIRSMAX, ALLPAIRSAVERAGE }; // Similarity measure.





class OntologyFactory {

public:

  /*!
    This object cannot be created.
  */
  OntologyFactory() = delete;
  ~OntologyFactory() = delete;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Information Factory.
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  [[nodiscard]] static std::shared_ptr <const InformationInterface> createInformation( const std::shared_ptr<const TermAnnotation>& annotation_ptr,
                                                                                       const std::shared_ptr<const GoGraphImpl>& graph_ptr,
                                                                                       InformationContentType information_type);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Similarity Factory.
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  [[nodiscard]] static std::shared_ptr<const SimilarityInterface> createSimilarity( const std::shared_ptr<const TermAnnotation>& annotation_ptr,
                                                                                    const std::shared_ptr<const GoGraphImpl>& graph_ptr,
                                                                                    SimilarityType similarity_type,
                                                                                    InformationContentType information_type = InformationContentType::CONTENT);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Similarity Factory.
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  [[nodiscard]] static std::shared_ptr<const SetSimilarityInterface> createSetSimilarity( const std::shared_ptr<const TermAnnotation>& annotation_ptr,
                                                                                          const std::shared_ptr<const GoGraphImpl>& graph_ptr,
                                                                                          SetSimilarityType set_similarity_type,
                                                                                          SimilarityType similarity_type = SimilarityType::LIN,
                                                                                          InformationContentType information_type = InformationContentType::CONTENT);

  // Call with the following; BESTMATCHAVERAGE, AVERAGEBESTMATCH, ALLPAIRSMAX, ALLPAIRSAVERAGE.
  [[nodiscard]] static std::shared_ptr<const SetSimilarityInterface> createSetSimilarity( const std::shared_ptr<const SimilarityInterface>& similarity_ptr,
                                                                                          SetSimilarityType set_similarity_type);

};



} // namespace


#endif //KOL_SIMILARITYFACTORYIMPL_H
