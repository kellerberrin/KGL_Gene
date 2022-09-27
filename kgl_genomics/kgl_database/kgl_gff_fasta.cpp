//
// Created by kellerberrin on 3/10/17.
//


#include "kel_exec_env.h"
#include "kel_utility.h"
#include "kgl_gff_fasta.h"
#include "kgl_sequence_base.h"

#include <functional>
#include <charconv>


namespace kgl = kellerberrin::genome;



std::shared_ptr<kgl::GenomeReference> kgl::ParseGffFasta::readFastaFile( const std::string& organism, const std::string& fasta_file_name) {

  std::vector<ReadFastaSequence> fasta_sequences;
  if (not readFastaFile(fasta_file_name, fasta_sequences)) {

    ExecEnv::log().critical("Could not read genome fasta file: {}", fasta_file_name);

  }

  std::shared_ptr<GenomeReference> genome_db_ptr(std::make_shared<GenomeReference>(organism));
  for (auto const& sequence : fasta_sequences) {

    StringDNA5 DNA5sequence(sequence.fastaSequence()); // convert to alphabet DNA5.
    std::shared_ptr<DNA5SequenceContig> sequence_ptr(std::make_shared<DNA5SequenceContig>(std::move(DNA5sequence)));
    const std::string contig_id = sequence.fastaId();

    if (not genome_db_ptr->addContigSequence(contig_id, sequence.fastaDescription(), sequence_ptr)) {

      ExecEnv::log().error("addContigSequence(), Attempted to add duplicate contig; {}", contig_id);

    }

    ExecEnv::log().info("Fasta Contig id: {}; Sequence length: {}; Description: {}", contig_id, sequence_ptr->length(), sequence.fastaDescription());

  }

  return genome_db_ptr;

}


bool kgl::ParseGffFasta::writeFastaFile( const std::string& fasta_file_name,
                                         const std::vector<WriteFastaSequence>& fasta_sequences) {

  return true;

}

bool kgl::ParseGffFasta::readFastaFile( const std::string& fasta_file_name,
                                        std::vector<ReadFastaSequence>& fasta_sequences) {

  // Simple state parser tokens.
  enum class ParserToken { FIND_FASTA_ID, FIND_FASTA_ID_NO_READ, PROCESS_FASTA_DATA, FASTA_PARSE_ERROR, READ_IO_EOF };

  // Open input file. Plain text or compressed.

  std::optional<std::unique_ptr<BaseStreamIO>> fasta_stream_opt = BaseStreamIO::getReaderStream(fasta_file_name);
  if (fasta_stream_opt) {

    ExecEnv::log().info("ParseGffFasta::readFastaFile; Opened fasta file: {} for processing", fasta_file_name);

  } else {

    ExecEnv::log().error("ParseGffFasta::readFastaFile; I/O error; could not open fasta file: {}", fasta_file_name);
    return false;

  }

  try {

    std::pair<std::string, std::string> fasta_id_comment;
    std::vector<std::unique_ptr<std::string>> fasta_lines;
    // Initial parser state.
    ParserToken parser_token = ParserToken::FIND_FASTA_ID;
    std::unique_ptr<std::string> record_ptr;

    do {

      if (parser_token != ParserToken::FIND_FASTA_ID_NO_READ) {

        IOLineRecord line_record = fasta_stream_opt.value()->readLine();
        if (line_record) {

          record_ptr = std::move(line_record.value().second);

        } else {

          parser_token = ParserToken::READ_IO_EOF;

        }

      }

      switch (parser_token) {

        case ParserToken::FIND_FASTA_ID_NO_READ:
        case ParserToken::FIND_FASTA_ID: {

        // Skip line if zero-sized, whitespace or comment.

          parser_token = ParserToken::FIND_FASTA_ID;
          if (record_ptr->empty() or std::isspace(record_ptr->front()) != 0 or record_ptr->front() == FASTA_COMMENT_) {

          break;

        }

        // Expect to find '>' char, complain and signal an error if not found.
        if (record_ptr->front() != FASTA_ID_) {

          ExecEnv::log().error("ParseGffFasta::readFastaFile; Expected to find fasta id line, found: {}", *record_ptr);
          parser_token = ParserToken::FASTA_PARSE_ERROR;
          break;

        }

        // Remove the first char '>'
        std::string id_line = *record_ptr;
        id_line.erase(id_line.begin());
        // Split into ID and comment on whitespace.
        fasta_id_comment = Utility::firstSplit(id_line);
        fasta_id_comment.second = Utility::trimEndWhiteSpace(fasta_id_comment.second);
        if (fasta_id_comment.first.empty()) {

          ExecEnv::log().error("Fasta file: {} has zero sized fasta id", fasta_file_name);
          parser_token = ParserToken::FASTA_PARSE_ERROR;
          break;

        }

        // Remove any fasta lines from the fasta data vector.
        fasta_lines.clear();
        parser_token = ParserToken::PROCESS_FASTA_DATA;

      }
      break;

      case ParserToken::PROCESS_FASTA_DATA: {

        //  If zero-sized or whitespace or comment or '>' then the fasta data block is complete.
        // Store the resultant festa record.
        if (record_ptr->empty() or std::isspace(record_ptr->front()) != 0 or record_ptr->front() == FASTA_COMMENT_ or record_ptr->front() == FASTA_ID_) {

          // Suppress IO, we already have an ID line.
          parser_token = ParserToken::FIND_FASTA_ID_NO_READ;
          // Create a fasta record.
          fasta_sequences.push_back(createFastaSequence(fasta_id_comment.first, fasta_id_comment.second, fasta_lines));

        } else {

          // Move the fasta data line into a vector of fasta lines.
          fasta_lines.push_back(std::move(record_ptr));

        }

      }
      break;

      case ParserToken::FASTA_PARSE_ERROR: {

        // Look for the next fasta id block.
        if (record_ptr) {

          if (record_ptr->front() == FASTA_ID_) {

            parser_token = ParserToken::FIND_FASTA_ID;

          }

        }

      }
      break;

      case ParserToken::READ_IO_EOF: {

        // Create a fasta record if data available.
        if (not fasta_id_comment.first.empty() and not fasta_lines.empty()) {

          fasta_sequences.push_back(createFastaSequence(fasta_id_comment.first, fasta_id_comment.second, fasta_lines));

        }

      }
      break;

      } // switch

    } while (parser_token != ParserToken::READ_IO_EOF);

  } catch(...) {

    ExecEnv::log().error("ParseGffFasta::readFastaFile; Unexpected IO exception");
    return false;

  }

  return true;

}


kgl::ReadFastaSequence kgl::ParseGffFasta::createFastaSequence( const std::string& fasta_id,
                                                                const std::string& fasta_comment,
                                                                const std::vector<std::unique_ptr<std::string>>& fasta_lines) {

  std::unique_ptr<std::string> fasta_data(std::make_unique<std::string>());
  size_t fasta_size{0};

  for (const auto& line : fasta_lines) {

    if (not line) {

      ExecEnv::log().error("ParseGffFasta::createFastaSequence; invalid fasta data line");

    } else {

      fasta_size += line->size();

    }

  }

  if (fasta_size == 0) {

    ExecEnv::log().error("ParseGffFasta::createFastaSequence; zero-sized fasta data");

  }

  fasta_data->reserve(fasta_size);

  for (const auto& line : fasta_lines) {

    if (not line) {

      ExecEnv::log().error("ParseGffFasta::createFastaSequence; invalid fasta data line");

    } else {

      std::copy(line->begin(), line->end(), std::back_inserter(*fasta_data));

    }

  }

  return {std::string(fasta_id), std::string(fasta_comment), std::move(fasta_data)};

}



void kgl::ParseGffFasta::readGffFile( const std::string &gff_file_name,
                                      kgl::GenomeReference& genome_db) {

  readGffFile(gff_file_name);

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void kgl::ParseGffFasta::readGffFile(const std::string& file_name) {

  FileDataIO file_io;
  std::vector<std::unique_ptr<GffRecord>> gff_records;
  size_t record_counter{0};

  if (not file_io.commenceIO(file_name)) {

    ExecEnv::log().critical("ParseGffFasta::readGffFile; I/O error; could not open file: {}", file_io.fileName());

  }

  while (true) {

    auto line_record = file_io.readIORecord();
    if (not line_record) break;

    auto const& [line_count, line_ptr] = line_record.value();

    const std::string& record_str = *line_ptr;

    if (not record_str.empty()) {

      if (record_str[0] == GFF_COMMENT_) {

        continue;  // Skip comment lines.

      }

    }

    auto [parse_result, gff_record_ptr] = parseGff3Record(std::move(line_record.value().second));
    if (not parse_result) {

      ExecEnv::log().error("ParseGffFasta::readGffFile; Bad row field format on line number: {}, text: {}", line_count, record_str);
      continue;

    }

    gff_records.push_back(std::move(gff_record_ptr));

    ++record_counter;

  }

  ExecEnv::log().info("ParseGffFasta::readGffFile; Parsed: {} lines, GFF3 records: {}", record_counter, gff_records.size());


}


std::pair<bool, std::unique_ptr<kgl::GffRecord>> kgl::ParseGffFasta::parseGff3Record(std::unique_ptr<std::string>&& gff_line) {

  std::unique_ptr<GffRecord> gff_record_ptr(std::make_unique<GffRecord>());
  std::vector<std::string_view> row_fields = Utility::view_tokenizer(*gff_line, GFF3_FIELD_DELIM_);

  if (row_fields.size() != GFF3_FIELD_COUNT_) {

    ExecEnv::log().error("ParseGffFasta::parseGff3Record; Bad field count: {}, text: {}", row_fields.size(), *gff_line);
    return {false, std::move(gff_record_ptr)};

  }

  bool parse_result{true};

  const std::string_view& id_field = row_fields[GFF3_ID_FIELD_IDX_];
  parse_result = parse_result and gff_record_ptr->id(id_field);

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
  const std::vector<std::string_view> tag_items = Utility::view_tokenizer(tag_field, GFF3_TAG_FIELD_DELIMITER_);

  std::vector<std::pair<std::string_view, std::string_view>> tag_value_vec;
  for (auto const& item : tag_items) {

    std::vector<std::string_view> tag_name = Utility::view_tokenizer(item, GFF3_TAG_ITEM_DELIMITER_);
    if (tag_name.size() != GFF3_ITEM_TAG_NAME_) {

      ExecEnv::log().error("ParseGffFasta::parseGff3Record; Bad 'tag=name' sub field: {}", item);
      parse_result = false;

    } else {

      tag_value_vec.emplace_back(std::pair(tag_name[0], tag_name[1]));

    }

  }

  parse_result = parse_result and gff_record_ptr->tagValues(tag_value_vec);

  return {parse_result, std::move(gff_record_ptr)};

}



std::shared_ptr<kgl::GenomeReference> kgl::ParseGffFasta::readFastaGffFile( const std::string& organism,
                                                                            const std::string& fasta_file_name,
                                                                            const std::string& gff_file_name ) {

  std::shared_ptr<GenomeReference> genome_db_ptr = readFastaFile(organism, fasta_file_name);
  readGffFile(gff_file_name, *genome_db_ptr);
  return genome_db_ptr;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// GffRecord
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::GffRecord::id(const std::string_view& id_txt) {

  if (id_txt.empty() or id_txt == MISSING_VALUE) {

    ExecEnv::log().error("GffRecord::id, feature record missing id, id text: {}", id_txt);
    id_.clear();
    return false;

  }

  id_ = id_txt;
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

  id_ = type_txt;
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

    phase_ = NO_PHASE;
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


bool kgl::GffRecord::tagValues(const std::vector<std::pair<std::string_view, std::string_view>>& tag_value_pairs) {

  bool result{true};
  for (auto const& [tag, value] : tag_value_pairs) {

    auto [iter, insert_result] = tagValues_.try_emplace(std::string(tag), std::string(value));
    if (not insert_result) {

      ExecEnv::log().error("GffRecord::tagValues; duplicate feature [tag, value] pair: [{}, {}]", tag, value);
      result = false;

    }

  }

  return result;

}
