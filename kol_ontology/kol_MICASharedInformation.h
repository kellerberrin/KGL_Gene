/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/

#ifndef KGL_MICA_SHARED_INFORMATION
#define KGL_MICA_SHARED_INFORMATION

#include "kol_InformationInterface.h"
#include "kol_GoGraph.h"

namespace kellerberrin::ontology {


/*! \class MICASharedInformation
	\brief A class to calculate shared information as the most informative common ancestor (MICA)

	This class calculates shared information using the most informative common ancestor (MICA).
	 The MICA is a term that is also known as the minimum subsumer.

	 This shared information method forms the basis of 3 information content measures
	 put forward by Lord el al.

    P. W. Lord, R. D. Stevens, A. Brass, and C. A. Goble, 
	 "Semantic similarity measures as tools for exploring the gene ontology,"
	 Pac Symp Biocomput, pp. 601-12, 2003.

*/
class MICASharedInformation : public InformationInterface {

public:
  //! A constructor
  /*!
    Creates the MICASharedInformation class
  */
  MICASharedInformation( const std::shared_ptr<const GoGraph> &graph_ptr,
                         const std::shared_ptr<const InformationInterface> &ic_map_ptr)
      : graph_ptr_(graph_ptr), ic_map_ptr_(ic_map_ptr) {}

  ~MICASharedInformation() override = default;

  //! A method for calculating the shared infromation between two concepts.
  /*!
    This method returns the shared information between two concepts.
  */
  [[nodiscard]] double sharedInformation(const std::string &termA, const std::string &termB) const override;

  //! An interface method for returning the shared information of a single terms,or information content
  /*!
    This method privdes a mechanism for returing a term's information content.
  */
  [[nodiscard]] double termInformation(const std::string &term) const override;

  //! An interface method for returning the maximum information content for a term
  /*!
    This method provides the absolute max information content within a corpus for normalization purposes.
  */
  [[nodiscard]] double maxInformationContent(const std::string &term) const override;

  //! Method to test if the ids exist and have the same ontology in the map
  /*!
    Return true the ids are found and the same ontology, false if not
  */
  [[nodiscard]] bool validateTerms(const std::string &id_termA, const std::string &id_termB) const override {

    return ic_map_ptr_->validateTerms(id_termA, id_termB);

  }

private:

  std::shared_ptr<const GoGraph> graph_ptr_;
  std::shared_ptr<const InformationInterface> ic_map_ptr_;

};

} // namespace


#endif
