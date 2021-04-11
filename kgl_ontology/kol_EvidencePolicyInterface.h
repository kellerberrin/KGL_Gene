/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_EVIDENCE_POLICY_INTERFACE
#define KGL_EVIDENCE_POLICY_INTERFACE

#include "kol_GoEnums.h"

namespace kellerberrin::ontology {

/*! \class EvidencePolicyInterface
	\brief An interface to check evidence codes for GO annotations

	This is interface is used to create parsers which will only use a specific set
	  of evidence codes when parsing annotations
*/
class EvidencePolicyInterface {

public:

  EvidencePolicyInterface() = default;

  virtual ~EvidencePolicyInterface() = default;

  //! A pure virtual method to test if an evidence code is allowed
  /*!
    This pure virtual method requires any subclass to imlement an isAllowed
      method to enforce the evidence pollicy.
  */
  [[nodiscard]] virtual bool isAllowed(const GO::EvidenceCode &evidenceCode) const = 0;

  /*!
    This pure virtual method copies an evidence policy.
  */
  [[nodiscard]] virtual std::unique_ptr<const EvidencePolicyInterface> clone() const = 0;

  /*!
    An invalid policy contains the error code 'GO::EvidenceCode::ECODE_ERROR'
  */
  [[nodiscard]] virtual bool isValid() const = 0;


};

} // namespace

#endif
