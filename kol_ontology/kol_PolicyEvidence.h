/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KOL_ALLOWED_SET_EVIDENCE_POLICY
#define KOL_ALLOWED_SET_EVIDENCE_POLICY

#include <vector>
#include "kol_GoEnums.h"
#include "kol_OntologyTypes.h"


namespace kellerberrin::ontology {

/*! \class PolicyEvidence
	\brief A class to allow only a set of evidence codes for annotations

	A class to allow only certain evidence codes in the go graph. It uses a set of 
	enums to restric the types of evidence codes considered for annotations.

*/
class PolicyEvidence {

public:

  //! A constructor
  /*!
    Creates the default PolicyEvidence, subscribes to all valid evidence types.
  */
  PolicyEvidence() { addEvidenceSet(GO::getAllValidEvidence()); }
  ~PolicyEvidence() = default;

  //! A parameterized constructor
  /*!
    Creates the PolicyEvidence using a list(vector) of evidence codes to allow
  */
  explicit PolicyEvidence(const std::vector<GO::EvidenceCode> &evidenceCodes) { addEvidenceSet(evidenceCodes); }


  //! a method to test if an evidence code is allowed or not
  /*!
    tests if the evidence is allowed. Overridden to fulfill the PolicyEvidenceInterface
  */
  [[nodiscard]] bool isAllowed(const GO::EvidenceCode &evidenceCode) const { return evidence_set_.contains(evidenceCode); }

  //! a method to add a evidence to the set of evidence codes allowed
  /*!
    adds a evidence to the set of evidence codes allowed by setting its mapped value to true
  */
  void addEvidenceSet(const std::vector<GO::EvidenceCode>& evidenceCodes)
  {

    for (auto const & evidence_code : evidenceCodes) {

      evidence_set_.insert(evidence_code);

    }

  }

  //! An invalid policy is empty or contains an invalid error code.
  /*!
    Determines if the Policy is valid
  */

  [[nodiscard]] bool isValid() const { return not evidence_set_.empty() and not evidence_set_.contains(GO::EvidenceCode::ECODE_ERROR); }


private:
  /*! The set of valid evidence codes */
  OntologySetType<GO::EvidenceCode> evidence_set_;

};


} // namespace

#endif
