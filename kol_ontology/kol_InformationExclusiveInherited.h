/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KOL_EXCLUSIVELY_INHERITED_SHARED_INFORMATION
#define KOL_EXCLUSIVELY_INHERITED_SHARED_INFORMATION

#include "kol_InformationContent.h"
#include "kol_GoGraph.h"



namespace kellerberrin::ontology {

/*! \class InformationExclusiveInherited
	\brief A class to calculate shared information in linear time after Zhang and Lai.

	Shu-Bo Zhang and Jian-Huang Lai. Semantic Similarity measurement between gene ontology
	  terms based on exclusively inherited shared information. Gene 558 (2015) 108-117.

*/
class InformationExclusiveInherited : public InformationInterface {

public:

  //! A constructor
  /*!
    Creates the InformationCoutoGraSM class
  */
  InformationExclusiveInherited(const std::shared_ptr<const GoGraphImpl> &graph_ptr,
                                const std::shared_ptr<const InformationContent> &ic_map_ptr)
      : graph_ptr_(graph_ptr), ic_map_ptr_(ic_map_ptr) {}

  ~InformationExclusiveInherited() override = default;


  //! An method for returning the shared information of two terms
  /*!
    This method returns the mean information content of the frontier ancestors
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

  std::shared_ptr<const GoGraphImpl> graph_ptr_;
  std::shared_ptr<const InformationContent> ic_map_ptr_;

  //! A method for determining the common disjunctive ancestors
  /*!
    This method returns the common disjunctive ancestors for two terms
  */
  [[nodiscard]] OntologySetType<std::string> getCommonDisjointAncestors(const std::string &termC1,
                                                                        const std::string &termC2) const;


};

} // namespace

#endif
