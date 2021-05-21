//
// Created by kellerberrin on 20/5/21.
//

#ifndef KOL_TERMINFORMATIONCONTENTIMPL_H
#define KOL_TERMINFORMATIONCONTENTIMPL_H



#include "kol_GoGraph.h"
#include "kol_AnnotationData.h"
#include "kol_TermInformationInterface.h"


namespace kellerberrin::ontology {


/*! \class TermInformationContentMap
	\brief A class to calculate the information content of a GO term.

	This class provides a map that returns the information content of a GO term. This
	  class is used by Information Content methods.

*/


using TermProbOntMap = OntologyMapType<std::string, std::pair<double, GO::Ontology>>;

class TermInformationContentImpl : public TermInformationInterface {

public:

  //! Accessor for probablities vector
  /*!
    Get the vector of values
  */

  [[nodiscard]] const TermProbOntMap &getValues() const { return probability_map_; }

  //! Method to test if the id exists in the map
  /*!
    Return true the ids are found and the same ontology, false if not
  */
  [[nodiscard]] bool validateTerms(const std::string &id_termA, const std::string &id_termB) const override;

  //! Mapping function to return the value mapped by key
  /*!
  Get the value mapped by the given key. A specified function for the [] operator
  */

  [[nodiscard]] double termInformation(const std::string &term_id) const override;

  //-------------------------------------------------------------------

  //! Get the maximum information content for an ontology class
  /*!
    This function returns the the maximum information content for an ontology class
  */

  [[nodiscard]] double maxInformationContent(const std::string& term_id) const override;

  //! Public method for calculating the most informative common ancestor value
  /*!
    This method searches the sets to determine the most informative ancestor.
  */

  [[nodiscard]] double sharedInformation(const std::string& go_termA, const std::string& go_termB, const GoGraph &graph) const override;

protected:

  //! A default constructor
  /*!
    This constructor creates an empty IC map. Should not be used.
  */
  TermInformationContentImpl() = default;
  ~TermInformationContentImpl() override = default;


  TermProbOntMap probability_map_;
  const static constexpr double BAD_INFO_VALUE_{0.0};
  double max_bp_information{0.0};
  double max_mf_information{0.0};
  double max_cc_information{0.0};

  void convertProbtoIC();
  [[nodiscard]] double getEfficientMICA( const OntologySetType<std::string> &smaller_set,
                                         const OntologySetType<std::string> &larger_set) const;
  [[nodiscard]] double getMaxInformation(GO::Ontology ontology) const;
  [[nodiscard]] double getRootCount(const std::string& root_id) const;

};



} // namespace









#endif //KOL_TERMINFORMATIONCONTENTIMPL_H
