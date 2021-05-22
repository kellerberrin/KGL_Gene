/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_SHARED_INFORMATION_INTERFACE
#define KGL_SHARED_INFORMATION_INTERFACE

#include <string>


namespace kellerberrin::ontology {

/*! \class InformationInterface
	\brief An interface class to define shared information calculations

	This class defines the interface  for shared information calculations. Pure virtual methods
	  require that shared information methods implement these.

*/
class InformationInterface {

public:

  InformationInterface() = default;

  virtual ~InformationInterface() = default;

  //! A pure virtual method for returning the shared information of two terms
  /*!
    This pure virtual method requires any shared information class to implement this method.
  */
  [[nodiscard]] virtual double sharedInformation(const std::string &termA, const std::string &termB) const = 0;

  //! A pure virtual method for returning the shared information of a single terms,or information content
  /*!
    This pure virtual method privdes a mechanism for returing a term's information content.
  */
  [[nodiscard]] virtual double termInformation(const std::string &term) const = 0;

  //! A pure virtual method for returning the maximum information content for a term
  /*!
    This pure virtual method requires any shared information class to implement this method.
    This method provides the absolute max information content with in a corpus for normalization purposes.
  */
  [[nodiscard]] virtual double maxInformationContent(const std::string &term) const = 0;

  //! Method to test if the ids exist and have the same ontology in the map
  /*!
    Return true the ids are found and the same ontology, false if not
  */
  [[nodiscard]] virtual bool validateTerms(const std::string &termA, const std::string &termB) const = 0;

};


} // namespace

#endif
