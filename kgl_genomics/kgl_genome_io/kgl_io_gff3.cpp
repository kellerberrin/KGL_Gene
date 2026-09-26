//
// Created by kellerberrin on 28/09/22.
//

#include "kel_basic_io.h"
#include "kel_exec_env.h"
#include "kel_utility.h"
#include "kgl_io_gff3.h"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void kgl::ParseGff3::readGffFile( const std::string &gff_file_name, kgl::GenomeReference& genome_db) {

  std::map<std::string, size_t> type_count;

  bool result = parseGffFile(gff_file_name, [&](std::unique_ptr<GffRecord>&& record_ptr) {

    ++type_count[record_ptr->type()];

    if (not parseGffRecord(genome_db, *record_ptr)) {

      ExecEnv::log().warn("ParseGff3::readGffFile; Error parsing feature in Contig: {}", record_ptr->contig());

    }

  });

  if (not result) {

    ExecEnv::log().warn("ParseGff3::readGffFile; One or more GFF3 records could not be parsed in file: {}", gff_file_name);

  }

  // Generate some feature statistics.
  for (auto const& [type, count] : type_count) {

    ExecEnv::log().info("ParseGff3::readGffFile; Feature Type: {}, Count: {}", type, count);

  }

}


std::pair<bool, std::vector<std::unique_ptr<kgl::GffRecord>>> kgl::ParseGff3::readGffFile(const std::string& file_name) {

  std::vector<std::unique_ptr<GffRecord>> gff_records;

  bool result = parseGffFile(file_name, [&](std::unique_ptr<GffRecord>&& record_ptr) {

    gff_records.push_back(std::move(record_ptr));

  });

  return {result, std::move(gff_records)};

}


bool kgl::ParseGff3::parseGffFile(const std::string& file_name,
                                  const std::function<void(std::unique_ptr<GffRecord>&&)>& record_sink) {

  bool result{true};
  size_t record_counter{0};

  std::optional<std::unique_ptr<BaseStreamIO>> gff_stream_opt = BaseStreamIO::getStreamIO(file_name);
  if (not gff_stream_opt) {

    ExecEnv::log().critical("ParseGff3::parseGffFile; I/O error; could not open file: {}", file_name);

  }

  ExecEnv::log().info("ParseGff3::parseGffFile; Opened GFF3 file: {} for processing", file_name);

  while (true) {

    // Get the line record.
    auto line_record = gff_stream_opt.value()->readLine();

    // Terminate on EOF
    if (line_record.EOFRecord()) break;

    // Get the line data.
    auto [line_count, record_str] = line_record.getLineData();

    // Strip a trailing carriage return so that CRLF files are treated as LF files.
    if (not record_str.empty() and record_str.back() == '\r') {

      record_str.pop_back();

    }

    // Check for empty string.
    if (record_str.empty()) {

      ExecEnv::log().warn("ParseGff3::parseGffFile; unexpected zero length line found at parser Line: {}", line_count);
      continue;

    }

    // Skip comments, but stop at the embedded FASTA directive (the remainder is sequence, not GFF).
    if (record_str.front() == GFF_COMMENT_) {

      if (record_str.starts_with(GFF3_FASTA_DIRECTIVE_)) {

        ExecEnv::log().info("ParseGff3::parseGffFile; '##FASTA' directive found on line: {}; embedded sequence ignored", line_count);
        break;

      }

      continue;  // Skip comment lines.

    }

    // Parse the gff3 line.
    auto [parse_result, gff_record_ptr] = parseGff3Record(record_str);

    if (not parse_result) {

      ExecEnv::log().error("ParseGff3::parseGffFile; Bad row field format on line number: {}, Line text: {}", line_count, record_str);
      result = false;
      continue;

    }

    record_sink(std::move(gff_record_ptr));

    ++record_counter;

  }

  ExecEnv::log().info("ParseGff3::parseGffFile; Parsed: {} GFF3 records", record_counter);

  return result;

}


std::pair<bool, std::unique_ptr<kgl::GffRecord>> kgl::ParseGff3::parseGff3Record(const std::string& gff_line) {

  std::unique_ptr<GffRecord> gff_record_ptr(std::make_unique<GffRecord>());
  std::vector<std::string_view> row_fields = Utility::viewTokenizer(gff_line, GFF3_FIELD_DELIM_);

  if (row_fields.size() != GFF3_FIELD_COUNT_) {

    ExecEnv::log().error("ParseGff3::parseGff3Record; Bad field count: {}, text: {}", row_fields.size(), gff_line);
    return {false, std::move(gff_record_ptr)};

  }

  bool parse_result{true};

  // Evaluate every field (no short-circuit) so that malformed lines report all their defects.
  const bool contig_result = gff_record_ptr->contig(row_fields[GFF3_CONTIG_FIELD_IDX_]);
  parse_result = parse_result and contig_result;

  const bool source_result = gff_record_ptr->source(row_fields[GFF3_SOURCE_FIELD_IDX_]);
  parse_result = parse_result and source_result;

  const bool type_result = gff_record_ptr->type(row_fields[GFF3_TYPE_FIELD_IDX_]);
  parse_result = parse_result and type_result;

  const bool start_result = gff_record_ptr->convertStartOffset(row_fields[GFF3_START_FIELD_IDX_]);
  parse_result = parse_result and start_result;

  const bool end_result = gff_record_ptr->convertEndOffset(row_fields[GFF3_END_FIELD_IDX_]);
  parse_result = parse_result and end_result;

  const bool score_result = gff_record_ptr->score(row_fields[GFF3_SCORE_FIELD_IDX_]);
  parse_result = parse_result and score_result;

  const bool strand_result = gff_record_ptr->strand(row_fields[GFF3_STRAND_FIELD_IDX_]);
  parse_result = parse_result and strand_result;

  const bool phase_result = gff_record_ptr->phase(row_fields[GFF3_PHASE_FIELD_IDX_]);
  parse_result = parse_result and phase_result;

  const std::string_view& tag_field = row_fields[GFF3_TAG_FIELD_IDX_];
  const std::vector<std::string_view> tag_items = Utility::viewTokenizer(tag_field, GFF3_TAG_FIELD_DELIMITER_);

  std::vector<std::pair<std::string_view, std::string_view>> tag_value_vec;
  for (auto const& item : tag_items) {

    // A missing attribute value (".") or a void item (trailing or doubled ';') is legal and skipped.
    if (item.empty() or item == GffRecord::MISSING_VALUE) {

      continue;

    }

    std::vector<std::string_view> tag_name = Utility::viewTokenizer(item, GFF3_TAG_ITEM_DELIMITER_);
    if (tag_name.size() != GFF3_TAG_ITEM_FIELD_COUNT_) {

      ExecEnv::log().error("ParseGff3::parseGff3Record; Bad 'tag=name' sub field: {}", item);
      parse_result = false;

    } else {

      tag_value_vec.emplace_back(tag_name[0], tag_name[1]);

    }

  }

  const bool attributes_result = gff_record_ptr->attributes(tag_value_vec);
  parse_result = parse_result and attributes_result;

  return {parse_result, std::move(gff_record_ptr)};

}


bool kgl::ParseGff3::parseGffRecord(GenomeReference& genome_db, const GffRecord& gff_record) {
  // Get the attributes.
  // Get (or construct) the feature ID.
  std::vector<kgl::FeatureIdent_t> feature_id_vec;
  kgl::FeatureIdent_t feature_id;
  if (not gff_record.attributes().getIds(feature_id_vec)) {

    // Construct an id
    feature_id = gff_record.type() + std::to_string(gff_record.begin());
    ExecEnv::log().warn("ParseGff3::parseGffRecord; 'ID' key not found; ID: {} generated", feature_id);

  } else if (feature_id_vec.size() > 1) {

    ExecEnv::log().warn("ParseGff3::parseGffRecord; Gff feature Has {} 'ID' values, choosing first value", feature_id_vec.size());

    feature_id = feature_id_vec.front();

  } else {

    feature_id = feature_id_vec.front();

  }
  // Get a pointer to the contig_ref_ptr.
  std::optional<std::shared_ptr<const kgl::ContigReference>> contig_opt = genome_db.getContigSequence(gff_record.contig());
  if (not contig_opt) {

    ExecEnv::log().error("ParseGff3::parseGffRecord; Could not find contig_ref_ptr: {}", gff_record.contig());
    return false;

  }

  // Check the record interval is well formed (GFF3 requires start <= end).
  if (gff_record.begin() > gff_record.end()) {

    ExecEnv::log().warn("ParseGff3::parseGffRecord; Feature: {} has invalid interval [{}, {})", feature_id, gff_record.begin(), gff_record.end());
    return false;

  }

  // Check that the type field "CDS" also has a valid phase.
  if (gff_record.type() == Feature::CDS_TYPE_ and gff_record.phase() == GffRecord::INVALID_PHASE) {

    ExecEnv::log().error("ParseGff3::parseGffRecord; Mis-match between valid phase and CDS record type");
    return false;

  }

  FeatureSequence sequence (gff_record.begin(), gff_record.end(), gff_record.strand(), gff_record.phase());
  std::shared_ptr<kgl::Feature> feature_ptr;

  // Create feature objects according to type.
  // Switch on hashed type strings for convenience.
  switch(Utility::hash64(gff_record.type())) {

    case Utility::hash64(GeneFeature::CODING_GENE_):
    case Utility::hash64(GeneFeature::PROTEIN_CODING_GENE_): // Alias GFF types for a Gene.
    case Utility::hash64(GeneFeature::NCRNA_GENE_):
      feature_ptr = std::make_shared<GeneFeature>(feature_id, gff_record.type(), contig_opt.value(), sequence);
      break;

    case Utility::hash64(Feature::CDS_TYPE_):
      feature_ptr = std::make_shared<Feature>(feature_id, gff_record.type(), Feature::CDS_TYPE_, contig_opt.value(), sequence);
      break;

    case Utility::hash64(Feature::MRNA_TYPE_):
      feature_ptr = std::make_shared<Feature>(feature_id, gff_record.type(), Feature::MRNA_TYPE_, contig_opt.value(), sequence);
      break;

    case Utility::hash64(Feature::UTR5_TYPE_):
      feature_ptr = std::make_shared<Feature>(feature_id, gff_record.type(), Feature::UTR5_TYPE_, contig_opt.value(), sequence);
      break;

    case Utility::hash64(Feature::UTR3_TYPE_):
      feature_ptr = std::make_shared<Feature>(feature_id, gff_record.type(), Feature::UTR3_TYPE_, contig_opt.value(), sequence);
      break;

    case Utility::hash64(Feature::TSS_TYPE_):
      feature_ptr = std::make_shared<Feature>(feature_id, gff_record.type(), Feature::TSS_TYPE_, contig_opt.value(), sequence);
      break;

    default:
      feature_ptr = std::make_shared<Feature>(feature_id, gff_record.type(), gff_record.type() , contig_opt.value(), sequence);
      break;

  }

  // Add in the attributes.
  feature_ptr->setAttributes(gff_record.attributes());
  // Annotate the contig_ref_ptr.
  std::shared_ptr<kgl::ContigReference> mutable_contig_ptr = std::const_pointer_cast<kgl::ContigReference>(contig_opt.value());
  bool result = mutable_contig_ptr->addContigFeature(feature_ptr);

  if (not result) {

    ExecEnv::log().error("ParseGff3::parseGffRecord; Could not add duplicate feature: {} to contig_ref_ptr: {}", feature_id, gff_record.contig());
    return false;

  }

  return true;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// GffRecord
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::GffRecord::contig(const std::string_view& contig_txt) {

  if (contig_txt.empty() or contig_txt == MISSING_VALUE) {

    ExecEnv::log().error("GffRecord::contig, feature record missing contig, contig text: {}", contig_txt);
    contig_.clear();
    return false;

  }

  contig_ = contig_txt;
  return true;

}



bool kgl::GffRecord::source(const std::string_view& source_txt) {

  if (source_txt == MISSING_VALUE) {

    source_.clear();

  } else {

    source_ = source_txt;

  }

  return true;

}

bool kgl::GffRecord::type(const std::string_view& type_txt) {

  if (type_txt.empty() or type_txt == MISSING_VALUE) {

    ExecEnv::log().error("GffRecord::type, feature record missing type, type text: {}", type_txt);
    type_.clear();
    return false;

  }

  type_ = type_txt;
  std::transform(type_.begin(), type_.end(), type_.begin(), [](unsigned char chr) {
    return static_cast<char>(std::toupper(chr));
  });

  return true;

}

bool kgl::GffRecord::convertStartOffset(const std::string_view& offset_txt) {

  auto [ptr, ec] = std::from_chars(offset_txt.data(), offset_txt.data() + offset_txt.size(), begin_position_);
  if (ec != ERRC_SUCCESS or begin_position_ == 0) {

    ExecEnv::log().error("GffRecord::convertStartOffset; bad feature start offset text: {}", offset_txt);
    begin_position_ = INVALID_OFFSET;
    return false;

  } else {

    --begin_position_;  // Start offset is ZERO based.

  }

  return true;

}

bool kgl::GffRecord::convertEndOffset(const std::string_view& offset_txt) {

  auto [ptr, ec] = std::from_chars(offset_txt.data(), offset_txt.data() + offset_txt.size(), end_position_);
  if (ec != ERRC_SUCCESS or end_position_ == 0) {

    ExecEnv::log().error("GffRecord::convertEndOffset; bad feature end offset text: {}", offset_txt);
    end_position_ = INVALID_OFFSET;
    return false;

  }

  return true;

}

bool kgl::GffRecord::score(const std::string_view& score_txt) {

  if (score_txt.empty() or score_txt == MISSING_VALUE) {

    score_ = NO_SCORE;
    return true;

  }

  auto [ptr, ec] = std::from_chars(score_txt.data(), score_txt.data() + score_txt.size(), score_);
  if (ec != ERRC_SUCCESS) {

    ExecEnv::log().error("GffRecord::score; bad feature score text: {}", score_txt);
    score_ = NO_SCORE;
    return false;

  }

  return true;

}


bool kgl::GffRecord::phase(const std::string_view& phase_txt) {

  if (phase_txt.empty() or phase_txt == MISSING_VALUE) {

    phase_ = NO_PHASE;
    return true;

  }

  auto [ptr, ec] = std::from_chars(phase_txt.data(), phase_txt.data() + phase_txt.size(), phase_);
  if (ec != ERRC_SUCCESS or phase_ > MAX_PHASE) {

    ExecEnv::log().error("GffRecord::phase; bad feature phase text: {}", phase_txt);
    phase_ = INVALID_PHASE;
    return false;

  }

  return true;

}

bool kgl::GffRecord::strand(const std::string_view& strand_txt) {

  if (strand_txt.empty() or strand_txt == MISSING_VALUE) {

    strand_ = StrandSense::FORWARD;
    return true;

  }

  if (strand_txt == STRAND_FORWARD_CHAR) {

    strand_ = StrandSense::FORWARD;
    return true;

  }

  if (strand_txt == STRAND_REVERSE_CHAR) {

    strand_ = StrandSense::REVERSE;
    return true;

  }

  ExecEnv::log().error("GffRecord::strand; Unexpected character: '{}' used to specify strand sense", strand_txt);
  strand_ = StrandSense::FORWARD;
  return false;

}


bool kgl::GffRecord::attributes(const std::vector<std::pair<std::string_view, std::string_view>>& tag_value_pairs) {

  for (auto const& [tag, value] : tag_value_pairs) {

    record_attributes_.insertAttribute(std::string(tag), std::string(value));

  }

  return true;

}
