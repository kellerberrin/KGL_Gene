
#ifndef KOL_PESQUITA_SIMGIC_SET_SIMILARITY
#define KOL_PESQUITA_SIMGIC_SET_SIMILARITY

#include "kol_SetSimilarityInterface.h"
#include "kol_InformationInterface.h"
#include "kol_GoGraph.h"


namespace kellerberrin::ontology {

/*! \class SetSimilarityPesquitaSimGIC
\brief A class to calculate Pesquita's SimGIC similarity between sets of go terms.

Pesquita, C., Faria, D., Bastos, H., Falcao, A., & Couto, F. (2007, July).
Evaluating GO-based semantic similarity measures. In Proc. 10th Annual
Bio-Ontologies Meeting (Vol. 37, No. 40, p. 38).

*/
class SetSimilarityPesquitaSimGIC : public SetSimilarityInterface {

public:
  //! Constructor
  /*!
  Creates the SetSimilarityPesquitaSimGIC class assigning the GoGraphImpl private member.
  */
  SetSimilarityPesquitaSimGIC(const std::shared_ptr<const GoGraph> &graph_ptr,
                              const std::shared_ptr<const InformationInterface> &ic_map_ptr)
      : graph_ptr_(graph_ptr), ic_map_ptr_(ic_map_ptr) {}

  ~SetSimilarityPesquitaSimGIC() override = default;

  //! A method for calculating term set to term set similarity for GO terms;
  /*!
  This method returns the best match average similarity.
  */
  [[nodiscard]] double calculateSimilarity( const OntologySetType<std::string> &row_terms,
                                            const OntologySetType<std::string> &column_terms) const override;

private:

  //! Pointer to the GoGraphImpl object
  /*!
  A reference to GO graph to be used.
  */
  std::shared_ptr<const GoGraph> graph_ptr_;

  //! The information content map
  /*!
  An information content map
  */
  std::shared_ptr<const InformationInterface> ic_map_ptr_;

};

} // namespace

#endif
