//
// Created by kellerberrin on 3/8/21.
//

#include "kgl_bio_pmid_parser.h"
#include "kel_utility.h"



namespace kgl = kellerberrin::genome;

kgl::BioPMIDRecord::BioPMIDRecord(const BioPMIDRecord& record) {

  pmid_id = record.pmid_id;
  bio_type = record.bio_type;
  bio_id = record.bio_id;
  bio_text = record.bio_text;

}


// Efficent move/copy operations for the Bio PMID record.

kgl::BioPMIDRecord::BioPMIDRecord(  std::string&& pmid_id_rval,
                                    std::string&& bio_type_rval,
                                    std::string&& bio_id_rval,
                                    std::string&& bio_text_rval) noexcept {

  pmid_id = pmid_id_rval;
  bio_type = bio_type_rval;
  bio_id = bio_id_rval;
  bio_text = bio_text_rval;

}



kgl::BioPMIDRecord::BioPMIDRecord(BioPMIDRecord&& record) noexcept {

  pmid_id = std::move(record.pmid_id);
  bio_type = std::move(record.bio_type);
  bio_id = std::move(record.bio_id);
  bio_text = std::move(record.bio_text);

}

kgl::BioPMIDRecord& kgl::BioPMIDRecord::operator=(BioPMIDRecord&& record) noexcept {

  pmid_id = std::move(record.pmid_id);
  bio_type = std::move(record.bio_type);
  bio_id = std::move(record.bio_id);
  bio_text = std::move(record.bio_text);

  return *this;

}



std::vector<std::string> kgl::BioPMIDResource::entrezPMID(const std::string& entrez_id) const {

  std::vector<std::string> pmid_vector;

  auto [lower_bound, upper_bound] = entrez_pmid_map_.equal_range(entrez_id);

  while (lower_bound != upper_bound) {

    auto const& [entrez_key, pmid_record] = *lower_bound;

    pmid_vector.push_back(pmid_record.pmid_id);

    ++lower_bound;

  }

  return pmid_vector;

}


std::vector<std::string> kgl::BioPMIDResource::diseaseMeSHPMID(const std::string& disease_mesh_id) const {

  std::vector<std::string> pmid_vector;

  auto [lower_bound, upper_bound] = disease_pmid_map_.equal_range(disease_mesh_id);

  while (lower_bound != upper_bound) {

    auto const& [mesh_key, pmid_record] = *lower_bound;

    pmid_vector.push_back(pmid_record.pmid_id);

    ++lower_bound;

  }

  return pmid_vector;

}


bool kgl::ParseBioPMID::parseBioPMIDFile(const std::string& file_name) {


  auto parsed_record_ptr = parseFlatFile(file_name, SquareTextParser::DELIMITER_TSV);

  if (parsed_record_ptr->getRowVector().size() < MINIMUM_ROW_COUNT_) {

    ExecEnv::log().error("ParseBioPMID::parseBioPMIDFile; Row count: {} for file: {} is below minimum",
                         parsed_record_ptr->getRowVector().size(), file_name);
    return false;

  }

  if (not parsed_record_ptr->checkRowSize(COLUMN_COUNT_)) {

    ExecEnv::log().error("ParseBioPMID::parseBioPMIDFile; Not all rows have expected column count: {} for file: {}",
                         COLUMN_COUNT_, file_name);
    return false;

  }

  ExecEnv::log().info("Begin Parsing BIO PMID Resource for file: {}", file_name);

  for (auto&& row_vector :  parsed_record_ptr->getRowVector()) {

    std::string bio_type = std::move(row_vector[BIO_TYPE_OFFSET_]);
    if (bio_type == DISEASE_TAG) {

      std::string bio_id = std::move(row_vector[BIO_ID_OFFSET_]);
      if (bio_id.empty()) {

        ExecEnv::log().warn("ParseBioPMID::parseBioPMIDFile; Bio ID key empty");
        continue;

      }

      std::string bio_key = bio_id;
      BioPMIDRecord bio_pmid_record(std::move(row_vector[PMID_OFFSET_]), std::move(bio_type), std::move(bio_id), std::move(row_vector[BIO_TEXT_OFFSET_]));
      disease_pmid_map_.emplace(std::move(bio_key), std::move(bio_pmid_record));

    } else if (bio_type == ENTREZ_GENE_TAG) {

      std::string bio_id = std::move(row_vector[BIO_ID_OFFSET_]);
      if (bio_id.empty()) {

        ExecEnv::log().warn("ParseBioPMID::parseBioPMIDFile; Bio ID key empty");
        continue;

      }

      std::string bio_key = bio_id;
      BioPMIDRecord bio_pmid_record(std::move(row_vector[PMID_OFFSET_]), std::move(bio_type), std::move(bio_id), std::move(row_vector[BIO_TEXT_OFFSET_]));
      entrez_pmid_map_.emplace(std::move(bio_key), std::move(bio_pmid_record));

    }

  }

  ExecEnv::log().info("ParseBioPMID::parseBioPMIDFile; Parsed: {} [PMID, Disease], Parsed: {} [PMID, Entrez_Gene] records from file: {}",
                      disease_pmid_map_.size() , entrez_pmid_map_.size(), file_name);

  return true;

}
