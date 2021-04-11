#ifndef KGL_RESNIK_SIMILARITY
#define KGL_RESNIK_SIMILARITY

#include <cmath>
#include <string>

#include "TermSimilarityInterface.h"
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
  ResnikSimilarity(const std::shared_ptr<const GoGraph> &goGraph, const std::shared_ptr<const TermInformationContentMap> &icMap)
      : _goGraph(goGraph), _icMap(icMap) {}

  ~ResnikSimilarity() override = default;

  //! A method for calculating term-to-term similarity for GO terms using Resnik similarity
  /*!
    This method returns the Resnik similarity or the information content of the most informative common ancestor.
  */
  [[nodiscard]] double calculateTermSimilarity(const std::string &goTermA, const std::string &goTermB) const override {
    //if the terms do not exit return 0.0 similarity
    if (not _icMap->hasTerm(goTermA) or not _icMap->hasTerm(goTermB)) {

      return 0.0;

    }

    //if not from same ontology, return 0;
    if (_goGraph->getTermOntology(goTermA) != _goGraph->getTermOntology(goTermB)) {

      return 0.0;

    }

    //Create 2 sets of term + ancestors
    OntologySetType<std::string> ancestorsA = _goGraph->getSelfAncestorTerms(goTermA);
    OntologySetType<std::string> ancestorsB = _goGraph->getSelfAncestorTerms(goTermB);

    return _icMap->getMICAinfo(ancestorsA, ancestorsB);

  }

  //! A method for calculating term-to-term similarity for GO terms using Normalized Resnik similarity
  /*!
    This method returns the Resnik similarity divided by the maximum possible similarity
  */
  [[nodiscard]] double calculateNormalizedTermSimilarity(const std::string &goTermA, const std::string &goTermB) const override {
    //call base similarity
    double resnik = calculateTermSimilarity(goTermA, goTermB);

    if (resnik <= 0) {

      return 0.0;

    }

    //select the correct ontology normalization factor
    GO::Ontology ontology = _goGraph->getTermOntology(goTermA);
    double ontology_info;
    switch (ontology) {

      case GO::Ontology::BIOLOGICAL_PROCESS:
        ontology_info = _icMap->getMinBP();
        break;

      case GO::Ontology::MOLECULAR_FUNCTION:
        ontology_info = _icMap->getMinMF();
        break;

      case GO::Ontology::CELLULAR_COMPONENT:
        ontology_info = _icMap->getMinCC();
        break;

      default:
      case GO::Ontology::ONTO_ERROR:
        ontology_info = 0.0;
        break;
    }

    if (ontology_info <= 0.0) {

      return 0.0;

    }

    double maxIC = -1.0 * std::log(ontology_info);
    return resnik / maxIC;

  }

private:


  std::shared_ptr<const GoGraph> _goGraph;
  std::shared_ptr<const TermInformationContentMap> _icMap;

};



} // namespace


#endif
