//
// Created by kellerberrin on 23/5/21.
//

#include "kol_GoEnums.h"


namespace kol = kellerberrin::ontology;


//! A method for returning the ontology code based on string
/*!
  This method takes a string and returns the proper enum
*/
kol::GO::Ontology kol::GO::ontologyStringToCode(const std::string &code) {

  for (auto const&[text, ontology] : ontology_text) {
    //return the index if matching, cast as enum
    if (code == text) {

      return ontology;

    }

  }
  //return error enum if not found
  return Ontology::ONTO_ERROR;

}

//! A method for returning a human readable string from the ontology code
/*!
  This method takes an ontology enum value and returns a string
*/
std::string kol::GO::ontologyToString(Ontology onto) {

  for (auto const&[text, ontology] : ontology_text) {

    if (ontology == onto) {

      return text;

    }

  }

  return ONTOLOGY_ERROR_TEXT;

}

// Retrieve evidence codes by evidence type.
const std::vector<kol::GO::EvidenceCode> kol::GO::getEvidenceType(EvidenceType evidence_type) {

  std::vector<EvidenceCode> evidence_codes;
  for (auto const& [text, code, type] : evidence_text) {

    if (type == evidence_type) {

      evidence_codes.push_back(code);

    }

  }

  return evidence_codes;

}

// Retrieve all valid evidence codes.
const std::vector<kol::GO::EvidenceCode> kol::GO::getAllValidEvidence() {

  std::vector<EvidenceCode> evidence_codes;
  for (auto const& [text, code, type] : evidence_text) {

    if (code != EvidenceCode::NR and code != EvidenceCode::ECODE_ERROR) {

      evidence_codes.push_back(code);

    }

  }

  return evidence_codes;

}

//! A method for converting evidence code strings to enums
/*!
  This method takes a string representing the evidence code and converts it to an enum.
*/
kol::GO::EvidenceCode kol::GO::evidenceStringToCode(const std::string &text_code) {

  for (auto const&[text, code, type] : evidence_text) {
    //return the index if matching, cast as enum
    if (text_code == text) {

      return code;

    }

  }
  //return evidence error if not found
  return EvidenceCode::ECODE_ERROR;

}

//! A method for returning a human readable string from an evidence code
/*!
  This method takes an evidence code enum value and returns a string
*/
std::string kol::GO::evidenceToString(EvidenceCode evidence) {

  for (auto const&[text, code, type] : evidence_text) {
    //return the index if matching, cast as enum
    if (evidence == code) {

      return text;

    }

  }
  //return evidence error if not found
  return EVIDENCE_ERROR_TEXT;

}



//! A method to convert relationship codes from string to enum
/*!
  This method converts the string representation of a relationship to an enum.
*/
kol::GO::Relationship kol::GO::relationshipStringToCode(const std::string &code) {

  for (auto const&[text, relation] : relationship_text) {
    //return the index if matching, cast as enum
    if (code == text) {

      return relation;

    }

  }
  //return error code if not found
  return Relationship::REL_ERROR;

}

//! A method for returning a human readable string from a Relationship
/*!
  This method takes an evidence code enum pointer value and returns a string
*/
std::string kol::GO::relationshipToString(Relationship relationship) {

  for (auto const&[text, relation] : relationship_text) {
    //return the index if matching, cast as enum
    if (relation == relationship) {

      return text;

    }

  }
  //return error code if not found
  return RELATIONSHIP_ERROR_TEXT;

}
