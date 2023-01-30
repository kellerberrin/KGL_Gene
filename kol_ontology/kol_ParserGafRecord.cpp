//
// Created by kellerberrin on 31/5/21.
//

#include "kel_exec_env.h"
#include "kel_utility.h"
#include "kol_ParserGafRecord.h"


namespace kol = kellerberrin::ontology;



bool kol::GAFRecord::parseGafRecord(const std::string& gaf_line_record) {

  std::vector<std::string> field_vec = Utility::charTokenizer(gaf_line_record, FIELD_SEPARATOR_);

  if (field_vec.size() != GO_FIELD_COUNT_) {

    ExecEnv::log().error("GAFRecord::parseGafRecord; Unexpected, Gaf record should have: {} fields, parsed: {} fields", GO_FIELD_COUNT_, field_vec.size());
    ExecEnv::log().error("GAFRecord::parseGafRecord; Gaf record : {}", gaf_line_record);
    return false;

  } else {

    DB_ident_ = std::move(field_vec[FIELD_OFFSET_DB_ID_]);              // required (1)
    gene_uniprot_id_ = std::move(field_vec[FIELD_OFFSET_UNIPROT_ID_]);            // required (1) key
    gene_symbolic_id_ = std::move(field_vec[FIELD_OFFSET_GENE_ID_]);	      // required (1)
    qualifier_ = std::move(field_vec[FIELD_OFFSET_QUALIFIER_]);          // optional (n)
    go_term_id_ = std::move(field_vec[FIELD_OFFSET_GO_ID_]);       // required (n)
    DB_reference_ = std::move(field_vec[FIELD_OFFSET_DB_REFERENCE_]); // required (n)
    evidence_code_ = GO::evidenceStringToCode(field_vec[FIELD_OFFSET_EVIDENCE_CODE_]);	    // required (n)
    if (evidence_code_ == GO::EvidenceCode::ECODE_ERROR) {

      ExecEnv::log().error("GAFRecord::parseGafRecord; Unknown/Invalid evidence code: {}", field_vec[FIELD_OFFSET_EVIDENCE_CODE_]);
      ExecEnv::log().error("GAFRecord::parseGafRecord; Gaf record : {}", gaf_line_record);
      return false;

    }
    with_from_ = std::move(field_vec[FIELD_OFFSET_WITH_FROM_]);	        // optional (n)
    const std::string& ontology_code = std::move(field_vec[FIELD_OFFSET_ONTOLOGY_]);             // required (n) biological (P), cellular (C) or molecular (F)

    if (ontology_code == ONTOLOGY_BP_CODE) {
      ontology_ = GO::Ontology::BIOLOGICAL_PROCESS;
    } else if (ontology_code == ONTOLOGY_MF_CODE) {
      ontology_ = GO::Ontology::MOLECULAR_FUNCTION;
    } else if (ontology_code == ONTOLOGY_CC_CODE) {
      ontology_ = GO::Ontology::CELLULAR_COMPONENT;
    } else {
      ExecEnv::log().error("GAFRecord::parseGafRecord; bad ontology code: {}", ontology_code);
      ExecEnv::log().error("GAFRecord::parseGafRecord; Gaf record : {}", gaf_line_record);
      return false;
    }

    description_ = std::move(field_vec[FIELD_OFFSET_NAME_]);	      // optional (1)
    alt_symbolic_ref_ = std::move(field_vec[FIELD_OFFSET_SYNONYM_]);       // required (1)
    DB_object_type_ = std::move(field_vec[FIELD_OFFSET_TYPE_]);	  // required (n)
    taxon_ = std::move(field_vec[FIELD_OFFSET_TAXON_]);             // required (1)
    date_ = std::move(field_vec[FIELD_OFFSET_DATE_]);              // required (1)
    assigned_by_ = std::move(field_vec[FIELD_OFFSET_ASSIGNED_]);       // required (n)
    annotation_ = std::move(field_vec[FIELD_OFFSET_EXTENSION_]);        // optional (n)
    gene_product_ = std::move(field_vec[FIELD_OFFSET_PRODUCT_]);      // optional (n)

  }

  return true;

}


