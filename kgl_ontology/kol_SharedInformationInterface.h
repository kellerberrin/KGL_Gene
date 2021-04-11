/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_SHARED_INFORMATION_INTERFACE
#define KGL_SHARED_INFORMATION_INTERFACE

#include <string>


namespace kellerberrin::ontology {

/*! \class SharedInformationInterface
	\brief An interface class to define shared information calculations

	This class defines the interface  for shared information calculations. Pure virtual methods
	  require that shared information methods implement these.

*/
class SharedInformationInterface {

public:

  SharedInformationInterface() = default;

  virtual ~SharedInformationInterface() = default;

  //! A pure virtual method for returning the shared information of two terms
  /*!
    This pure virtual method requires any shared information class to implement this method.
  */
  [[nodiscard]] virtual double sharedInformation(const std::string &termA, const std::string &termB) const = 0;

  //! A pure virtual method for returning the shared information of a single terms,or information content
  /*!
    This pure virtual method privdes a mechanism for returing a term's infromation content.
  */
  [[nodiscard]] virtual double sharedInformation(const std::string &term) const = 0;

  //! A pure virtual method for returning the maximum information content for a term
  /*!
    This pure virtual method requires any shared information class to implement this method.
    This method provides the absolute max information content with in a corpus for normalization purposes.
  */
  [[nodiscard]] virtual double maxInformationContent(const std::string &term) const = 0;

  //! A pure virtual method for determining if a term can be found.
  /*!
    This pure virtual method requires any shared information class to implement this method.
    This method provides a method for client classes to determine if a term can be found by the method.
  */
  [[nodiscard]] virtual bool hasTerm(const std::string &term) const = 0;

  //! A pure virtual method for determining if the two terms are of like ontologies.
  /*!
    This pure virtual method requires any shared information class to implement this method.
    This method provides a method for client classes to determine if two terms are of the same ontology.
  */
  [[nodiscard]] virtual bool isSameOntology(const std::string &termA, const std::string &termB) const = 0;


};

} // namespace

#endif
