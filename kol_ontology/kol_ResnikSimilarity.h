#ifndef KGL_RESNIK_SIMILARITY
#define KGL_RESNIK_SIMILARITY

#include <cmath>
#include <string>

#include "kol_TermSimilarityInterface.h"
#include "kol_TermInformationContentMap.h"
#include "kol_GoGraph.h"


namespace kellerberrin::ontology {


/*! \class ResnikSimilarity
	\brief A class to calculate resnik similarity between 2 terms

	This class calculates Resnik similarity.
	Philip Resnik (1995). "Using information content to evaluate semantic similarity in a taxonomy". 
	  In Chris S. Mellish (Ed.). Proceedings of the 14th international joint conference on Artificial intelligence (IJCAI'95)

    P. W. Lord, R. D. Stevens, A. Brass, and C. A. Goble, 
	 "Semantic similarity measures as tools for exploring the gene ontology,"
	 Pac Symp Biocomput, pp. 601-12, 2003.

	maximun information content of all shared ancestors
	IC(MICA)

*/
class ResnikSimilarity : public TermSimilarityInterface {

public:

  //! A constructor
  /*!
    Creates the default(empty) StandardRelationshipPolicy
  */
  ResnikSimilarity( const std::shared_ptr<const GoGraph> &graph_ptr,
                    const std::shared_ptr<const TermInformationInterface> &ic_map_ptr)
      : graph_ptr_(graph_ptr), ic_map_ptr_(ic_map_ptr) {}

  ~ResnikSimilarity() override = default;

  //! A method for calculating term-to-term similarity for GO terms using Resnik similarity
  /*!
    This method returns the Resnik similarity or the information content of the most informative common ancestor.
  */
  [[nodiscard]] double calculateTermSimilarity(const std::string &goTermA, const std::string &goTermB) const override;

  //! A method for calculating term-to-term similarity for GO terms using Normalized Resnik similarity
  /*!
    This method returns the Resnik similarity divided by the maximum possible similarity
  */
  [[nodiscard]] double calculateNormalizedTermSimilarity(const std::string &goTermA, const std::string &goTermB) const override;

private:

  std::shared_ptr<const GoGraph> graph_ptr_;
  std::shared_ptr<const TermInformationInterface> ic_map_ptr_;

};



} // namespace


#endif
