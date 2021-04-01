/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef GO_ENUMS
#define GO_ENUMS

#include <string>


//! Number of evidence codes, 22 + 1 error code
static const constexpr size_t NUM_EVIDENCES = 23;

//! Number of relationships codes, 5 + 1 error code
static const constexpr size_t NUM_RELATIONSHIPS = 6;


//! GO namespaces
/*!
	This namespace is a set of static variables related to go terms and relationships.
*/
namespace GO{

	//! function that returns strings representing the root ontology term biological_process
  static const constexpr char* ROOT_TERM_BIOLOICAL_PROCESS = "GO:0008150";
	[[nodiscard]] std::string getRootTermBP(){ return ROOT_TERM_BIOLOICAL_PROCESS; }

	//! function that returns strings representing the root ontology term molecular_function
  static const constexpr char* ROOT_TERM_MOLECULAR_FUNCTION = "GO:0003674";
	[[nodiscard]] std::string getRootTermMF(){ return ROOT_TERM_MOLECULAR_FUNCTION; }
	
	//! function that returns strings representing the root ontology term cellular_component
  static const constexpr char* ROOT_TERM_CELLULAR_COMPONENT = "GO:0005575";
	[[nodiscard]] std::string getRootTermCC(){ return ROOT_TERM_CELLULAR_COMPONENT; }

	//! Ontology enum type
	/*!
		This enum defines a type for sub-ontologies.
	*/
	enum class Ontology {
		BIOLOGICAL_PROCESS=0,
		MOLECULAR_FUNCTION=1,
		CELLULAR_COMPONENT=2,
		ONTO_ERROR=3
	};

	//! Ontology enum strings
	/*!
		These are string representations of the sub-ontologies. Take from obo go files.
	*/

	static const constexpr char* ONTOLOGY_BIOLOGICAL_PROCESS_TEXT = "biological_process";
  static const constexpr char* ONTOLOGY_MOLECULAR_FUNCTION_TEXT = "molecular_function";
  static const constexpr char* ONTOLOGY_CELLULAR_COMPONENT_TEXT = "cellular_component";
  static const constexpr char* ONTOLOGY_ERROR_TEXT = "ONTOLOGY_ERROR";

  struct OntologyText { const char* text; Ontology onto; };
	static const constexpr OntologyText ontology_text[] = {

	    { ONTOLOGY_BIOLOGICAL_PROCESS_TEXT, Ontology::BIOLOGICAL_PROCESS},
      { ONTOLOGY_MOLECULAR_FUNCTION_TEXT, Ontology::MOLECULAR_FUNCTION},
      { ONTOLOGY_CELLULAR_COMPONENT_TEXT, Ontology::CELLULAR_COMPONENT},
      { ONTOLOGY_ERROR_TEXT, Ontology::ONTO_ERROR}

  };


	//! A method for returning the ontology code based on string
	/*!
		This method takes a string and returns the proper enum
	*/
	[[nodiscard]] Ontology ontologyStringToCode(const std::string& code) {

		for(auto const& [text, ontology] : ontology_text) {
			//return the index if matching, cast as enum
			if(code == text) {

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
[[nodiscard]] std::string ontologyToString(Ontology onto) {

    for(auto const& [text, ontology] : ontology_text) {

      if(ontology == onto) {

        return text;

      }

    }

    return ONTOLOGY_ERROR_TEXT;

	}

	//! Evcidence Code enum type
	/*!
		This enum defines a type for evidence codes.
		 Defined at http://www.geneontology.org/GO.evidence.shtml
	*/
	enum class EvidenceCode{
		//experimental
		EXP=0,
		IDA=1,
		IPI=2,
		IMP=3,
		IGI=4,
		IEP=5,

		//computationally assisted
		ISS=6,
		ISO=7,
		ISA=8,
		ISM=9,
		IGC=10,
		IBA=11,
		IBD=12,
		IKR=13,
		IRD=14,
		RCA=15,

		//author statement
		TAS=16,
		NAS=17,

		//Curator statement
		IC=18,
		ND=19,

		//automatically assigned
		IEA=20,

		//obsolete evidence code
		NR=21,
		ECODE_ERROR=22
	};

  enum class EvidenceType { EXPERIMENTAL, COMPUTATIONAL, AUTHOR, CURATOR, AUTO_ASSIGNED, OBSOLETE, ERROR };

	//! Relationship code enum strings
	/*!
		These are string representations of the evidence codes for an annotation.
	*/
  struct EvidenceText { const char* text; EvidenceCode code; EvidenceType type; };
  static const constexpr char* EVIDENCE_ERROR_TEXT = "EVIDENCE_CODE_ERROR";
	static const constexpr EvidenceText evidence_text[] = {

		//experimental
      { "EXP", EvidenceCode::EXP, EvidenceType::EXPERIMENTAL },
      { "IDA", EvidenceCode::IDA, EvidenceType::EXPERIMENTAL },
      { "IPI", EvidenceCode::IPI, EvidenceType::EXPERIMENTAL },
      { "IMP", EvidenceCode::IMP, EvidenceType::EXPERIMENTAL },
      { "IGI", EvidenceCode::IGI, EvidenceType::EXPERIMENTAL },
      { "IEP", EvidenceCode::IEP, EvidenceType::EXPERIMENTAL },

		//computationally assisted
      { "ISS", EvidenceCode::ISS, EvidenceType::COMPUTATIONAL },
      { "ISO", EvidenceCode::ISO, EvidenceType::COMPUTATIONAL },
      { "ISA", EvidenceCode::ISA, EvidenceType::COMPUTATIONAL },
      { "ISM", EvidenceCode::ISM, EvidenceType::COMPUTATIONAL },
      { "IGC", EvidenceCode::IGC, EvidenceType::COMPUTATIONAL },
      { "IBA", EvidenceCode::IBA, EvidenceType::COMPUTATIONAL },
      { "IBD", EvidenceCode::IBD, EvidenceType::COMPUTATIONAL },
      { "IKR", EvidenceCode::IKR, EvidenceType::COMPUTATIONAL },
      { "IRD", EvidenceCode::IRD, EvidenceType::COMPUTATIONAL },
      { "RCA", EvidenceCode::RCA, EvidenceType::COMPUTATIONAL },

		//author statement
      { "TAS", EvidenceCode::TAS, EvidenceType::AUTHOR },
      { "NAS", EvidenceCode::NAS, EvidenceType::AUTHOR },

		//Curator statement
      { "IC", EvidenceCode::IC, EvidenceType::CURATOR },
      { "ND", EvidenceCode::ND, EvidenceType::CURATOR },

		//automatically assigned
      { "IEA", EvidenceCode::IEA, EvidenceType::AUTO_ASSIGNED },

		//obsolete evidence code
      { "NR", EvidenceCode::NR, EvidenceType::OBSOLETE },

		//Error code
      { EVIDENCE_ERROR_TEXT, EvidenceCode::ECODE_ERROR, EvidenceType::ERROR }
	};

	//! A method for converting evidence code strings to enums
	/*!
		This method takes a string representing the evidence code and converts it to an enum.
	*/
	[[nodiscard]] EvidenceCode evidenceStringToCode(const std::string& text_code){

		for(auto const& [text, code, type] : evidence_text) {
			//return the index if matching, cast as enum
			if(text_code == text) {

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
	[[nodiscard]] std::string evidenceToString(EvidenceCode evidence){

    for(auto const& [text, code, type] : evidence_text) {
      //return the index if matching, cast as enum
      if (evidence == code) {

        return text;

      }

    }
    //return evidence error if not found
    return EVIDENCE_ERROR_TEXT;

  }


	//! Relationship codes enum
	/*!
		This enum represents the relationship codes for ontology edges.
	*/
	enum class Relationship {
		IS_A=0,
		PART_OF=1,
		REGULATES=2,
		POSITIVELY_REGULATES=3,
		NEGATIVELY_REGULATES=4,
		REL_ERROR=5
	};


//! Relationship code enum strings
/*!
  These strings represent the enum codes for each relationship.
*/

struct RelationshipText { const char* text; Relationship type; };
static const constexpr char* RELATIONSHIP_ERROR_TEXT = "RELATIONSHIP_ERROR";
static const constexpr RelationshipText relationship_text[] = {

    { "is_a", Relationship::IS_A },
    { "part_of", Relationship::PART_OF },
    { "regulates", Relationship::REGULATES },
    { "positively_regulates", Relationship::POSITIVELY_REGULATES },
    { "negatively_regulates", Relationship::NEGATIVELY_REGULATES },
    { RELATIONSHIP_ERROR_TEXT, Relationship::REL_ERROR }

	};

	//! A method to convert relationship codes from string to enum
	/*!
		This method converts the string representation of a relationship to an enum.
	*/
	[[nodiscard]] Relationship relationshipStringToCode(const std::string& code){

		for(auto const& [text, relation] : relationship_text) {
			//return the index if matching, cast as enum
			if(code == text){

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
	[[nodiscard]] std::string relationshipToString(Relationship relationship){

    for(auto const& [text, relation] : relationship_text) {
      //return the index if matching, cast as enum
      if(relation == relationship){

        return text;

      }

    }
    //return error code if not found
		return RELATIONSHIP_ERROR_TEXT;

	}

}
#endif
