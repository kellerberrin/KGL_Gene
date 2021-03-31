/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef EXPERIMENTAL_EVIDENCE_POLICY
#define EXPERIMENTAL_EVIDENCE_POLICY

#include <AllowedSetEvidencePolicy.hpp>

/*! \class ExperimentalEvidencePolicy
	\brief A class to allow experimental evidence codes for annotations

	A class to allow only experimental evidence codes for gene annotations. This class
	 extends the AllowedSetEvidencePolicy and adds the experimental evidence codes to
	 the allowed set.
*/
class ExperimentalEvidencePolicy: public AllowedSetEvidencePolicy{

public:
	//! A constructor
	/*!
		Creates the default(empty) AllowedSetEvidencePolicy
	*/
	ExperimentalEvidencePolicy() : AllowedSetEvidencePolicy() {

		addEvidence(GO::EvidenceCode::EXP);
		addEvidence(GO::EvidenceCode::IDA);
		addEvidence(GO::EvidenceCode::IPI);
		addEvidence(GO::EvidenceCode::IMP);
		addEvidence(GO::EvidenceCode::IGI);
		addEvidence(GO::EvidenceCode::IEP);

	}


  ~ExperimentalEvidencePolicy() override = default;

  [[nodiscard]] std::unique_ptr<const EvidencePolicyInterface> clone() const override { return std::make_unique<const ExperimentalEvidencePolicy>(); }

};
#endif
