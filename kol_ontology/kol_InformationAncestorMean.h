/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KOL_ANCESTOR_MEAN_SHARED_INFORMATION
#define KOL_ANCESTOR_MEAN_SHARED_INFORMATION

#include "kol_InformationContentImpl.h"
#include "kol_GoGraph.h"


namespace kellerberrin::ontology {

/*! \class InformationAncestorMean
	\brief A class to calculate shared infromation as the average information conent of all common ancestors

	This class calculates shared information by averaging the information content of all common ancestors.

	 This shared information method is used a baseline for comparison and may not be meaningful.

*/
class InformationAncestorMean : public InformationInterface {

public:
  //! A constructor
  /*!
    Creates the InformationAncestorMean class
  */
  InformationAncestorMean(const std::shared_ptr<const GoGraph> &graph_ptr,
                          const std::shared_ptr<const InformationContentImpl> &ic_map_ptr)
      : graph_ptr_(graph_ptr), ic_map_ptr_(ic_map_ptr) {}

  ~InformationAncestorMean() override = default;

  //! A method for calculating the shared infromation between two concepts.
  /*!
    This method returns the shared information between two concepts.
  */
  [[nodiscard]] double sharedInformation(const std::string &termA, const std::string &termB) const override;

  //! An interface method for returning the shared information of a single terms,or information content
  /*!
    This method privdes a mechanism for returing a term's infromation content.
  */
  [[nodiscard]] double termInformation(const std::string &term) const override { return ic_map_ptr_->termInformation(term); }


  //! An interface method for returning the maximum information content for a term
  /*!
    This method provides the absolute max information content within a corpus for normalization purposes.
  */
  [[nodiscard]] double maxInformationContent(const std::string &term) const override { return ic_map_ptr_->maxInformationContent(term); }

  //! Method to test if the ids exist and have the same ontology in the map
  /*!
    Return true the ids are found and the same ontology, false if not
  */
  [[nodiscard]] bool validateTerms(const std::string &id_termA, const std::string &id_termB) const override {

    return ic_map_ptr_->validateTerms(id_termA, id_termB);

  }


private:

  std::shared_ptr<const GoGraph> graph_ptr_;
  std::shared_ptr<const InformationContentImpl> ic_map_ptr_;


};

} // namespace

#endif
