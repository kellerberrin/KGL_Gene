//
// Created by kellerberrin on 3/8/21.
//

#include "kgl_bio_pmid_parser.h"
#include "kel_mt_buffer.h"
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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::set<std::string> kgl::BioPMIDMaps::entrezPMID(const std::string& entrez_id) const {

  std::set<std::string> pmid_set;

  auto [lower_bound, upper_bound] = entrez_pmid_map_.equal_range(entrez_id);

  while (lower_bound != upper_bound) {

    auto const& [entrez_key, pmid_record] = *lower_bound;

    pmid_set.insert(pmid_record.pmid_id);

    ++lower_bound;

  }

  return pmid_set;

}


std::set<std::string> kgl::BioPMIDMaps::selectDiseaseBioPMID( const std::vector<std::string>& or_MeSHlist,
                                                              const std::vector<std::string>& and_MeSHlist) const {

  std::set<std::string> or_pmid_set;
  for (auto const& mesh_term : or_MeSHlist) {

    auto pmid_set = diseaseMeSHPMID(mesh_term);
    for (auto const& pmid : pmid_set) {

      or_pmid_set.insert(pmid);

    }

  }

  if (and_MeSHlist.empty()) {

    return or_pmid_set;

  }

  std::set<std::string> and_pmid_set;
  for (auto const& mesh_term : and_MeSHlist) {

    auto pmid_set = diseaseMeSHPMID(mesh_term);
    for (auto const& pmid : pmid_set) {

      and_pmid_set.insert(pmid);

    }

  }

  std::set<std::string> result_set;
  for (const auto& pmid : or_pmid_set) {

    if (and_pmid_set.contains(pmid)) {

      result_set.insert(pmid);

    }

  }

  return result_set;

}


std::set<std::string> kgl::BioPMIDMaps::diseaseMeSHPMID(const std::string& disease_mesh_id) const {

  std::set<std::string> pmid_set;

  auto [lower_bound, upper_bound] = disease_pmid_map_.equal_range(disease_mesh_id);

  while (lower_bound != upper_bound) {

    auto const& [mesh_key, pmid_record] = *lower_bound;

    pmid_set.insert(pmid_record.pmid_id);

    ++lower_bound;

  }

  return pmid_set;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::ParseBioPMID::parseBioPMIDRecords(const std::string& file_name) {

  StreamMTBuffer file_io;
  size_t counter = 0;

  if (not file_io.open(file_name)) {

    ExecEnv::log().error("ParseBioPMID::parseBioPMIDRecords; I/O error; could not open file: {}", file_name);
    return false;

  }

  while (true) {

    auto line_record = file_io.readLine();
    if (line_record.EOFRecord()) {

      break; // EOF reached.

    }

    const std::string_view& record_view = line_record.getView();

    if (not record_view.empty()) {

      if (record_view[0] == COMMENT_) {

        continue;  // Skip comment lines.

      }

    }

    std::vector<std::string_view> row_fields = Utility::viewTokenizer(record_view, DELIMITER_);

    if (row_fields.size() != COLUMN_COUNT_) {

      ExecEnv::log().error("ParseBioPMID::parseBioPMIDRecords; Incorrect field count: {}, expected: {}, line: {}, text: {}",
                           row_fields.size(), COLUMN_COUNT_, line_record.lineCount(), record_view);
      continue;

    }

    if (not parseFields(row_fields)) {

      ExecEnv::log().error("ParseBioPMID::parseBioPMIDRecords; Problem parsing line: {}", line_record.lineCount());

    }

    ++counter;

    if (counter % REPORT_INTERVAL_ == 0) {

      ExecEnv::log().info("Parsed Bio concept Pubmed PMID Records: {}", counter);

    }

  }

  ExecEnv::log().info("Parsed: {} lines from text file: {}", counter, file_name);
  return true;

}


bool kgl::ParseBioPMID::parseFields(const std::vector<std::string_view>& field_views) {

  if (field_views[BIO_ID_OFFSET_].empty()) {

    ExecEnv::log().warn("ParseBioPMID::parseBioPMIDFile; Bio ID key empty");
    return false;

  }

  if (field_views[BIO_TYPE_OFFSET_] == DISEASE_TAG) {


    disease_pmid_map_.emplace(std::string(field_views[BIO_ID_OFFSET_]), BioPMIDRecord( std::string(field_views[PMID_OFFSET_]),
                                                                                       std::string(field_views[BIO_TYPE_OFFSET_]),
                                                                                       std::string(field_views[BIO_ID_OFFSET_]),
                                                                                       std::string(field_views[BIO_TEXT_OFFSET_])));

  } else if (field_views[BIO_TYPE_OFFSET_] == ENTREZ_GENE_TAG) {

    entrez_pmid_map_.emplace(std::string(field_views[BIO_ID_OFFSET_]), BioPMIDRecord( std::string(field_views[PMID_OFFSET_]),
                                                                                      std::string(field_views[BIO_TYPE_OFFSET_]),
                                                                                      std::string(field_views[BIO_ID_OFFSET_]),
                                                                                      std::string(field_views[BIO_TEXT_OFFSET_])));

  }

  return true;

}
