//
// Created by kellerberrin on 28/09/22.
//




#include "kel_exec_env.h"
#include "kel_utility.h"
#include "kgl_gff3.h"
#include "kel_mt_buffer.h"

#include <functional>
#include <charconv>


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void kgl::ParseGff3::readGffFile( const std::string &gff_file_name, kgl::GenomeReference& genome_db) {

  std::map<std::string, size_t> type_count;

  auto [result, record_ptr_vector] = readGffFile(gff_file_name);
  if (result) {

    for (auto& record_ptr : record_ptr_vector) {

      if (type_count.contains(record_ptr->type())) {

        ++type_count[record_ptr->type()];

      } else {

        type_count[record_ptr->type()] = 1;

      }

      if (not parseGffRecord(genome_db, *record_ptr)) {

        ExecEnv::log().warn("ParseGff3::readGffFile; Error parsing feature in Contig: {}", record_ptr->contig());

      }

    }

  }

  // Generate some feature statistics.
  for (auto const& [type, count] : type_count) {

    ExecEnv::log().info("ParseGffFasta::readGffFile; Feature Type: {}, Count: {}", type, count);

  }

}


std::pair<bool, std::vector<std::unique_ptr<kgl::GffRecord>>> kgl::ParseGff3::readGffFile(const std::string& file_name) {

  std::vector<std::unique_ptr<GffRecord>> gff_records;
  size_t record_counter{0};
  bool result{true};
  StreamMTBuffer file_io;

  if (not file_io.open(file_name)) {

    ExecEnv::log().critical("ParseGffFasta::readGffFile; I/O error; could not open file: {}", file_name);

  }

  ExecEnv::log().info("ParseGffFasta::readGffFile; Opened GFF3 file: {} for processing", file_name);

  while (true) {

    // Get the line record.
    auto line_record = file_io.readLine();

    // Terminate on EOF
    if (line_record.EOFRecord()) break;

    // Get the line data.
    auto const [line_count, record_str] = line_record.getLineData();

    // Skip comments
    if (record_str[0] == GFF_COMMENT_) {

      continue;  // Skip comment lines.

    }

    // Check for empty string
    if (record_str.empty()) {

      ExecEnv::log().warn("ParseGffFasta::readGffFile; unexpected zero length line found at parser Line: {}", line_count);
      continue;

    }

    // Parse the gff3 line.
    auto [parse_result, gff_record_ptr] = parseGff3Record(record_str);

    if (not parse_result) {

      ExecEnv::log().error("ParseGffFasta::readGffFile; Bad row field format on line number: {}, Line text: {}", line_count, record_str);
      continue;

    }

    gff_records.push_back(std::move(gff_record_ptr));

    result = result and parse_result;

    ++record_counter;

  }

  ExecEnv::log().info("ParseGffFasta::readGffFile; Parsed: {} lines, GFF3 records: {}", record_counter, gff_records.size());

  return {result, std::move(gff_records)};

}


std::pair<bool, std::unique_ptr<kgl::GffRecord>> kgl::ParseGff3::parseGff3Record(const std::string& gff_line) {

  std::unique_ptr<GffRecord> gff_record_ptr(std::make_unique<GffRecord>());
  std::vector<std::string_view> row_fields = Utility::viewTokenizer(gff_line, GFF3_FIELD_DELIM_);

  if (row_fields.size() != GFF3_FIELD_COUNT_) {

    ExecEnv::log().error("ParseGffFasta::parseGff3Record; Bad field count: {}, text: {}", row_fields.size(), gff_line);
    return {false, std::move(gff_record_ptr)};

  }

  bool parse_result{true};

  const std::string_view& contig_field = row_fields[GFF3_CONTIG_FIELD_IDX_];
  parse_result = parse_result and gff_record_ptr->contig(contig_field);

  const std::string_view& source_field = row_fields[GFF3_SOURCE_FIELD_IDX_];
  parse_result = parse_result and gff_record_ptr->source(source_field);

  const std::string_view& type_field = row_fields[GFF3_TYPE_FIELD_IDX_];
   parse_result = parse_result and gff_record_ptr->type(type_field);

  const std::string_view& start_field = row_fields[GFF3_START_FIELD_IDX_];
  parse_result = parse_result and gff_record_ptr->convertStartOffset(start_field);

  const std::string_view& end_field = row_fields[GFF3_END_FIELD_IDX_];
  parse_result = parse_result and gff_record_ptr->convertEndOffset(end_field);

  const std::string_view& score_field = row_fields[GFF3_SCORE_FIELD_IDX_];
  parse_result = parse_result and gff_record_ptr->score(score_field);

  const std::string_view& strand_field = row_fields[GFF3_STRAND_FIELD_IDX_];
  parse_result = parse_result and gff_record_ptr->strand(strand_field);

  const std::string_view& phase_field = row_fields[GFF3_PHASE_FIELD_IDX_];
  parse_result = parse_result and gff_record_ptr->phase(phase_field);

  const std::string_view& tag_field = row_fields[GFF3_TAG_FIELD_IDX_];
  const std::vector<std::string_view> tag_items = Utility::viewTokenizer(tag_field, GFF3_TAG_FIELD_DELIMITER_);

  std::vector<std::pair<std::string_view, std::string_view>> tag_value_vec;
  for (auto const& item : tag_items) {

    std::vector<std::string_view> tag_name = Utility::viewTokenizer(item, GFF3_TAG_ITEM_DELIMITER_);
    if (tag_name.size() != GFF3_ITEM_TAG_NAME_) {

      ExecEnv::log().error("ParseGffFasta::parseGff3Record; Bad 'tag=name' sub field: {}", item);
      parse_result = false;

    } else {

      tag_value_vec.emplace_back(std::pair(tag_name[0], tag_name[1]));

    }

  }

  parse_result = parse_result and gff_record_ptr->attributes(tag_value_vec);

  return {parse_result, std::move(gff_record_ptr)};

}


// Valgrind indicates memory leaks ocurring in this function.
// This could be due to the use of the FeatureSinkPtr function pointer to process
// the genomic features.
// todo: This memory leak needs further investigation.
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
  // Get a pointer to the contig.
  std::optional<std::shared_ptr<const kgl::ContigReference>> contig_opt = genome_db.getContigSequence(gff_record.contig());
  if (not contig_opt) {

    ExecEnv::log().error("ParseGff3::parseGffRecord; Could not find contig: {}", gff_record.contig());
    return false;

  }

  // Check that the type field "CDS" also has a valid phase.
  if (gff_record.type() == Feature::CDS_TYPE_ and gff_record.phase() == GffRecord::INVALID_PHASE) {

    ExecEnv::log().error("ParseGff3::parseGffRecord; Mis-match between valid phase and CDS record type");
    return false;

  }

  FeatureSequence sequence (gff_record.begin(), gff_record.end(), gff_record.strand());
  std::shared_ptr<kgl::Feature> feature_ptr;

  // Create feature objects according to type.
  // Switch on hashed type strings for convenience.
  switch(Utility::hash(gff_record.type())) {

    case Utility::hash(Feature::GENE_TYPE_):
    case Utility::hash(PROTEIN_CODING_GENE_): // Alias GFF types for a Gene.
    case Utility::hash(NCRNA_GENE_):
      feature_ptr = std::make_shared<GeneFeature>(feature_id, contig_opt.value(), sequence);
      break;

    case Utility::hash(Feature::CDS_TYPE_):
      feature_ptr = std::make_shared<CDSFeature>(feature_id, gff_record.phase(), contig_opt.value(), sequence);
      break;

    case Utility::hash(Feature::MRNA_TYPE_):
      feature_ptr = std::make_shared<Feature>(feature_id, Feature::MRNA_TYPE_, contig_opt.value(), sequence);
      break;

    case Utility::hash(Feature::UTR5_TYPE_):
      feature_ptr = std::make_shared<Feature>(feature_id, Feature::UTR5_TYPE_, contig_opt.value(), sequence);
      break;

    case Utility::hash(Feature::UTR3_TYPE_):
      feature_ptr = std::make_shared<Feature>(feature_id, Feature::UTR3_TYPE_, contig_opt.value(), sequence);
      break;

    case Utility::hash(Feature::TSS_TYPE_):
      feature_ptr = std::make_shared<Feature>(feature_id, Feature::TSS_TYPE_, contig_opt.value(), sequence);
      break;

    default:
      feature_ptr = std::make_shared<Feature>(feature_id, gff_record.type() , contig_opt.value(), sequence);
      break;

  }

  // Add in the attributes.
  feature_ptr->setAttributes(gff_record.attributes());
  // Annotate the contig.
  std::shared_ptr<kgl::ContigReference> mutable_contig_ptr = std::const_pointer_cast<kgl::ContigReference>(contig_opt.value());
  bool result = mutable_contig_ptr->addContigFeature(feature_ptr);

  if (not result) {

    ExecEnv::log().error("ParseGff3::parseGffRecord; Could not add duplicate feature: {} to contig: {}", feature_id, gff_record.contig());
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

    ExecEnv::log().error("GffRecord::id, feature record missing id, id text: {}", contig_txt);
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

    ExecEnv::log().error("GffRecord:type, feature record missing type, type text: {}", type_txt);
    type_.clear();
    return false;

  }

  type_ = type_txt;
  std::transform(type_.begin(), type_.end(), type_.begin(), ::toupper);

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

    ExecEnv::log().error("GffRecord::convertStartOffset; bad feature start offset text: {}", offset_txt);
    end_position_ = INVALID_OFFSET;
    return false;

  }

  return true;

}

bool kgl::GffRecord::score(const std::string_view& offset_txt) {

  if (offset_txt.empty() or offset_txt == MISSING_VALUE) {

    score_ = NO_SCORE;
    return true;

  }

  auto [ptr, ec] = std::from_chars(offset_txt.data(), offset_txt.data() + offset_txt.size(), score_);
  if (ec != ERRC_SUCCESS) {

    ExecEnv::log().error("GffRecord::score; bad feature score text: {}", offset_txt);
    score_ = NO_SCORE;
    return false;

  }

  return true;

}


bool kgl::GffRecord::phase(const std::string_view& offset_txt) {

  if (offset_txt.empty() or offset_txt == MISSING_VALUE) {

    phase_ = INVALID_PHASE;
    return true;

  }

  auto [ptr, ec] = std::from_chars(offset_txt.data(), offset_txt.data() + offset_txt.size(), phase_);
  if (ec != ERRC_SUCCESS or phase_ > MAX_PHASE) {

    ExecEnv::log().error("GffRecord::phase; bad feature phase text: {}", offset_txt);
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

  bool result{true};
  for (auto const& [tag, value] : tag_value_pairs) {

    record_attributes_.insertAttribute(std::string(tag), std::string(value));

  }

  return result;

}
