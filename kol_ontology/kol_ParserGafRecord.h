//
// Created by kellerberrin on 31/5/21.
//

#ifndef KOL_PARSER_GAFRECORD_H
#define KOL_PARSER_GAFRECORD_H


#include "kol_GoEnums.h"

#include <vector>
#include <string>

namespace kellerberrin::ontology {


class GAFRecord {

public:

  GAFRecord() = default;

  GAFRecord(const GAFRecord &) = default;

  ~GAFRecord() = default;

  // Getters
  [[nodiscard]] const std::string &DBIdent() const { return DB_ident_; }              // required
  [[nodiscard]] const std::string &geneUniprotId() const { return gene_uniprot_id_; }    // required
  [[nodiscard]] const std::string &geneSymbolicId() const { return gene_symbolic_id_; }  // required
  [[nodiscard]] const std::string &qualifier() const { return qualifier_; }        // optional
  [[nodiscard]] const std::string &goIdent() const { return go_term_id_; }          // required
  [[nodiscard]] const std::string &DBReference() const { return DB_reference_; }      // required
  [[nodiscard]] std::string evidenceCodeText() const { return GO::evidenceToString(evidence_code_); }      // required
  [[nodiscard]] GO::EvidenceCode evidenceCode() const { return evidence_code_; }      // required
  [[nodiscard]] const std::string &withFrom() const { return with_from_; }          // optional
  [[nodiscard]] GO::Ontology ontology() const { return ontology_; }             // required
  [[nodiscard]] std::string ontologyText() const { return GO::ontologyToString(ontology_); }             // required
  [[nodiscard]] const std::string &description() const { return description_; }        // optional
  [[nodiscard]] const std::string &altSymbolicRef() const { return alt_symbolic_ref_; } // required
  [[nodiscard]] const std::string &DBObjectType() const { return DB_object_type_; }      // required
  [[nodiscard]] const std::string &taxon() const { return taxon_; }              // required
  [[nodiscard]] const std::string &date() const { return date_; }               // required
  [[nodiscard]] const std::string &assignedBy() const { return assigned_by_; }       // required
  [[nodiscard]] const std::string &annotation() const { return annotation_; }         // optional
  [[nodiscard]] const std::string &geneProduct() const { return gene_product_; }       // optional

  // Setter.
  [[nodiscard]] bool parseGafRecord(const std::string &gaf_line_record);

private:

  constexpr static const size_t GO_FIELD_COUNT_ = 17;
  constexpr static const char FIELD_SEPARATOR_ = '\t';
  const static constexpr size_t FIELD_OFFSET_DB_ID_ = 0;
  const static constexpr size_t FIELD_OFFSET_UNIPROT_ID_ = 1;
  const static constexpr size_t FIELD_OFFSET_GENE_ID_ = 2;
  const static constexpr size_t FIELD_OFFSET_QUALIFIER_ = 3;
  const static constexpr size_t FIELD_OFFSET_GO_ID_ = 4;
  const static constexpr size_t FIELD_OFFSET_DB_REFERENCE_ = 5;
  const static constexpr size_t FIELD_OFFSET_EVIDENCE_CODE_ = 6;
  const static constexpr size_t FIELD_OFFSET_WITH_FROM_ = 7;
  const static constexpr size_t FIELD_OFFSET_ONTOLOGY_ = 8;
  const static constexpr char *ONTOLOGY_BP_CODE = "P";
  const static constexpr char *ONTOLOGY_MF_CODE = "F";
  const static constexpr char *ONTOLOGY_CC_CODE = "C";
  const static constexpr size_t FIELD_OFFSET_NAME_ = 9;
  const static constexpr size_t FIELD_OFFSET_SYNONYM_ = 10;
  const static constexpr size_t FIELD_OFFSET_TYPE_ = 11;
  const static constexpr size_t FIELD_OFFSET_TAXON_ = 12;
  const static constexpr size_t FIELD_OFFSET_DATE_ = 13;
  const static constexpr size_t FIELD_OFFSET_ASSIGNED_ = 14;
  const static constexpr size_t FIELD_OFFSET_EXTENSION_ = 15;
  const static constexpr size_t FIELD_OFFSET_PRODUCT_ = 16;


  std::string DB_ident_;              // required
  std::string gene_uniprot_id_;    // required
  std::string gene_symbolic_id_;    // required
  std::string qualifier_;          // optional
  std::string go_term_id_;          // required
  std::string DB_reference_;       // required
  GO::EvidenceCode evidence_code_;      // required
  std::string with_from_;          // optional
  GO::Ontology ontology_;             // required
  std::string description_;        // optional
  std::string alt_symbolic_ref_;  // required
  std::string DB_object_type_;      // required
  std::string taxon_;              // required
  std::string date_;               // required
  std::string assigned_by_;        // required
  std::string annotation_;         // optional
  std::string gene_product_;       // optional

};


} // namespace

#endif //KGL_KOL_PARSERGAFRECORD_H
