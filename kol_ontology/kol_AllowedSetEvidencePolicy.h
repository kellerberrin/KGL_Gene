/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_ALLOWED_SET_EVIDENCE_POLICY
#define KGL_ALLOWED_SET_EVIDENCE_POLICY

#include <vector>
#include "kol_GoEnums.h"
#include "kol_EvidencePolicyInterface.h"


namespace kellerberrin::ontology {

/*! \class AllowedSetEvidencePolicy
	\brief A class to allow only a set of evidence codes for annotations

	A class to allow only certain evidence codes in the go graph. It uses a set of 
	enums to restric the types of evidence codes considered for annotations.

*/
class AllowedSetEvidencePolicy : public EvidencePolicyInterface {

public:

  //! A constructor
  /*!
    Creates the default(empty) AllowedSetEvidencePolicy
  */
  AllowedSetEvidencePolicy() = default;

  AllowedSetEvidencePolicy(const AllowedSetEvidencePolicy &copy) = default;

  ~AllowedSetEvidencePolicy() override = default;

  //! A parameterized constructor
  /*!
    Creats the AllowedSetEvidencePolicy using a list(vector) of evidence codes to allow
  */
  explicit AllowedSetEvidencePolicy(const std::vector<GO::EvidenceCode> &evidenceCodes) {

    for (auto const &evidence :  evidenceCodes) {

      _evidenceMap[evidence] = true;

    }

  }


  //! A parameterized constructor
  /*!
    Clone create the DisallowedSetEvidencePolicy.
  */

  [[nodiscard]] std::unique_ptr<const EvidencePolicyInterface> clone() const override { return std::make_unique<const AllowedSetEvidencePolicy>(*this); }


  //! a method to test if an eviddence code is allowed or not
  /*!
    tests if the evidence is allowed. Overridden to fulfill the EvidencePolicyInterface
  */
  [[nodiscard]] bool isAllowed(const GO::EvidenceCode &evidenceCode) const override { return _evidenceMap.find(evidenceCode) != _evidenceMap.end(); }

  //! a method to add a evidence to the set of evidence codes allowed
  /*!
    adds a evidence to the set of evidence codes allowed by setting its mapped value to true
  */
  void addEvidence(GO::EvidenceCode evidenceCode) { _evidenceMap[evidenceCode] = true; }

  //! a method to add a evidence to the set of evidence codes allowed
  /*!
    adds a evidence to the set of evidence codes allowed by setting its mapped value to true
  */
  void addEvidence(const std::string &stringCode) {

    GO::EvidenceCode evidenceCode = GO::evidenceStringToCode(stringCode);

    if (evidenceCode != GO::EvidenceCode::ECODE_ERROR) {

      _evidenceMap[evidenceCode] = true;

    }

  }


  //! An invalid policy contains the error code 'GO::EvidenceCode::ECODE_ERROR'
  /*!
    Determines if the Policy is valid
  */

  [[nodiscard]] bool isValid() const override { return _evidenceMap.find(GO::EvidenceCode::ECODE_ERROR) == _evidenceMap.end(); }


private:
  //! a map of evidence codes to bool
  /*! maps an evidence code to a bool. Boost unordered map give constant time find. */
  OntologyMapType<GO::EvidenceCode, bool> _evidenceMap;

};


} // namespace

#endif
