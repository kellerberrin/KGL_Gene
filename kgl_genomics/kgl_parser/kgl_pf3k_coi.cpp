//
// Created by kellerberrin on 9/1/21.
//

#include "kgl_pf3k_coi.h"


namespace kgl = kellerberrin::genome;


kgl::Pf3kCOIResource::Pf3kCOIResource(std::string identifier, const SquareTextRows& square_text)
: ResourceBase(ResourceProperties::PF3K_COI_RESOURCE_ID_, std::move(identifier)) {

  if (not parseFlatFile(square_text)) {

    ExecEnv::log().error("Unable to parse input file for resource: {}", ResourceProperties::PF3K_COI_RESOURCE_ID_);

  }

  if (indexedFile().getHeaderMap().size() != FIELD_COUNT) {

    ExecEnv::log().error("Resource: {}, expected fields: {}, actual fields: {}",
    ResourceProperties::PF3K_COI_RESOURCE_ID_, Pf3kCOIResource::FIELD_COUNT, indexedFile().getHeaderMap().size());

  }

}


std::optional<size_t> kgl::Pf3kCOIResource::genomeCOI(const GenomeId_t& genome) const {

  const auto& genome_row = indexed_file_.getRow(genome);

  if (genome_row.empty()) {

    return std::nullopt;

  }

  if (genome_row.size() != FIELD_COUNT) {

    ExecEnv::log().warn("COI resource record for genome: {} contains fields: {}, expected fields: {}", genome, genome_row.size(), FIELD_COUNT);
    return std::nullopt;

  }

  return std::stoll(genome_row[COI_FIELD]);

}


bool kgl::Pf3kCOIParser::parseCOIPf3k(const std::string& file_name) {

  parsed_text_ptr_ = flat_file_parser_.parseFlatFile(file_name, SquareTextParser::DELIMITER_TSV);

  return true;

}
