/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef DISALLOWED_SET_EVIDENCE_POLICY
#define DISALLOWED_SET_EVIDENCE_POLICY

#include <vector>
#include <unordered_map>
#include <EvidencePolicyInterface.h>

/*! \class DisallowedSetEvidencePolicy
	\brief A class to allow only a set of evidence codes for annotations

	A class to allow only certain evidence codes in the go graph. It uses a set of 
	enums to restric the types of evidence codes considered for annotations.

*/
class DisallowedSetEvidencePolicy: public EvidencePolicyInterface{

public:
	
	//! A constructor
	/*!
		Creates the default(empty) DisallowedSetEvidencePolicy
	*/
	DisallowedSetEvidencePolicy() = default;
  DisallowedSetEvidencePolicy(const DisallowedSetEvidencePolicy& copy) = default;
  ~DisallowedSetEvidencePolicy() override = default;

	//! A parameterized constructor
	/*!
		Clone create the DisallowedSetEvidencePolicy.
	*/

  [[nodiscard]] std::unique_ptr<const EvidencePolicyInterface> clone() const override { return std::make_unique<const DisallowedSetEvidencePolicy>(*this); }

  //! a method to test if an evidence code is allowed or not
	/*!
		tests if the evidence is allowed. Overridden to fulfill the EvidencePolicyInterface.
		This class disallows a set and so returns false if the evidence code is found.
	*/
	[[nodiscard]] bool isAllowed(const GO::EvidenceCode& evidenceCode) const override { return _evidenceMap.find(evidenceCode) == _evidenceMap.end(); }

	//! a method to add a evidence to the set of evidence codes not allowed
	/*!
		adds a evidence to the set of evidence codes not allowed
	*/
	void addEvidence(const GO::EvidenceCode& evidenceCode) { _evidenceMap[evidenceCode] = true; }

	//! a method to add a evidence to the set of evidence codes allowed
	/*!
		adds a evidence to the set of evidence codes allowed by setting its mapped value to true
	*/
	void addEvidence(const std::string &stringCode){

		GO::EvidenceCode evidenceCode = GO::evidenceStringToCode(stringCode);

		if (evidenceCode != GO::EvidenceCode::ECODE_ERROR){

			_evidenceMap[evidenceCode] = true;

		}

	}


  //! An invalid policy contains the error code 'GO::EvidenceCode::ECODE_ERROR'
  /*!
    Determines if the Policy is valid
  */

  [[nodiscard]] bool isValid() const override { return _evidenceMap.find(GO::EvidenceCode::ECODE_ERROR) == _evidenceMap.end(); }


private:
	//! a map of evidence codes to boo
    /*! maps an evidence code to a bool. Boost unordered map give constant time find. */
	std::unordered_map<GO::EvidenceCode,bool> _evidenceMap;

};
#endif
