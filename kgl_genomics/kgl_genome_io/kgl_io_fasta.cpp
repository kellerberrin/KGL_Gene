//
// Created by kellerberrin on 28/09/22.
//

#include "kel_basic_io.h"
#include "kel_exec_env.h"
#include "kel_utility.h"
#include "kgl_io_fasta.h"

#include <cctype>
#include <fstream>
#include <string>
#include <utility>
#include <vector>



namespace kgl = kellerberrin::genome;


std::shared_ptr<kgl::GenomeReference> kgl::ParseFasta::readFastaFile( const std::string& organism, const std::string& fasta_file_name) {

  std::vector<ReadFastaSequence> fasta_sequences;
  if (not readFastaFile(fasta_file_name, fasta_sequences)) {

    ExecEnv::log().critical("ParseFasta::readFastaFile; Could not read genome fasta file: {}", fasta_file_name);

  }

  std::shared_ptr<GenomeReference> genome_db_ptr(std::make_shared<GenomeReference>(organism));
  for (auto const& sequence : fasta_sequences) {

    StringDNA5 DNA5sequence(sequence.fastaSequence()); // convert to alphabet DNA5.
    std::shared_ptr<DNA5SequenceLinear> sequence_ptr(std::make_shared<DNA5SequenceLinear>(std::move(DNA5sequence)));
    const std::string contig_id = sequence.fastaId();

    if (not genome_db_ptr->addContigSequence(contig_id, sequence.fastaDescription(), sequence_ptr)) {

      ExecEnv::log().error("ParseFasta::readFastaFile; addContigSequence(), Attempted to add duplicate contig_ref_ptr; {}", contig_id);

    }

    ExecEnv::log().info("ParseFasta::readFastaFile; Fasta Contig id: {}; Sequence length: {}; Description: {}", contig_id, sequence_ptr->length(), sequence.fastaDescription());

  }

  return genome_db_ptr;

}


bool kgl::ParseFasta::writeFastaFile( const std::string& fasta_file_name,
                                      const std::vector<WriteFastaSequence>& fasta_sequences) {

  std::ofstream fasta_file(fasta_file_name);

  if (not fasta_file.good()) {

    ExecEnv::log().error("ParseFasta::writeFastaFile; could not open file: {} for fasta file output", fasta_file_name);
    return false;

  } else {

    ExecEnv::log().info("ParseFasta::writeFastaFile; writing fasta record: {} to file: {}", fasta_sequences.size(), fasta_file_name);

  }

  for (auto const& fasta_record : fasta_sequences) {

    fasta_file << FASTA_ID_
               << fasta_record.fastaId()
               << " "
               << fasta_record.fastaDescription();

    auto nucleotide_view = fasta_record.fastaSequence().getStringView();

    for (size_t i = 0; i < nucleotide_view.size(); ++i) {

      if (i % FASTA_LINE_LENGTH_ == 0) {

        fasta_file << "\n";

      }

      fasta_file << nucleotide_view[i];

    }

    // Terminate the record. Without this the next '>id' header is appended to the last sequence line.
    fasta_file << "\n";

  }

  return fasta_file.good();

}

bool kgl::ParseFasta::isSkippableLine(std::string_view line) {

  return line.empty()
         or std::isspace(static_cast<unsigned char>(line.front())) != 0
         or line.front() == FASTA_COMMENT_;

}

bool kgl::ParseFasta::isFastaIdLine(std::string_view line) {

  return not line.empty() and line.front() == FASTA_ID_;

}

bool kgl::ParseFasta::readFastaFile( const std::string& fasta_file_name,
                                     std::vector<ReadFastaSequence>& fasta_sequences) {

  // Simple state parser tokens.
  enum class ParserToken { FIND_FASTA_ID, FIND_FASTA_ID_NO_READ, PROCESS_FASTA_DATA, FASTA_PARSE_ERROR, READ_IO_EOF };

  // Open input file. Plain text or compressed.

  std::optional<std::unique_ptr<BaseStreamIO>> fasta_stream_opt = BaseStreamIO::getStreamIO(fasta_file_name);
  if (fasta_stream_opt) {

    ExecEnv::log().info("ParseFasta::readFastaFile; Opened fasta file: {} for processing", fasta_file_name);

  } else {

    ExecEnv::log().error("ParseFasta::readFastaFile; I/O error; could not open fasta file: {}", fasta_file_name);
    return false;

  }

  try {

    std::pair<std::string, std::string> fasta_id_comment;
    std::vector<std::string> fasta_lines;
    // Initial parser state.
    ParserToken parser_token = ParserToken::FIND_FASTA_ID;
    std::pair<size_t, std::string> line_data;

    do {

      if (parser_token != ParserToken::FIND_FASTA_ID_NO_READ) {

        IOLineRecord line_record = fasta_stream_opt.value()->readLine();
        if (line_record.EOFRecord()) {

          parser_token = ParserToken::READ_IO_EOF;

        } else {

          line_data = line_record.getLineData();

        }

      }

      auto& [line_count, record_text] = line_data;

      // Strip a trailing carriage return so that CRLF files are treated as LF files.
      if (not record_text.empty() and record_text.back() == '\r') {

        record_text.pop_back();

      }

      switch (parser_token) {

        case ParserToken::FIND_FASTA_ID_NO_READ:
        case ParserToken::FIND_FASTA_ID: {

          // Skip line if zero-sized, whitespace or comment.

          parser_token = ParserToken::FIND_FASTA_ID;
          if (isSkippableLine(record_text)) {

            break;

          }

          // Expect to find '>' char, complain and signal an error if not found.
          if (not isFastaIdLine(record_text)) {

            ExecEnv::log().error("ParseFasta::readFastaFile; Expected to find fasta id line, found: {}", record_text);
            parser_token = ParserToken::FASTA_PARSE_ERROR;
            break;

          }

          // Remove the first char '>'
          std::string id_line = record_text;
          id_line.erase(id_line.begin());
          // Split into ID and comment on whitespace.
          fasta_id_comment = Utility::firstSplitChar(id_line);
          fasta_id_comment.second = Utility::trimEndWhiteSpace(fasta_id_comment.second);
          if (fasta_id_comment.first.empty()) {

            ExecEnv::log().error("ParseFasta::readFastaFile; Fasta file: {} has zero sized fasta id", fasta_file_name);
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
          // Store the resultant fasta record.
          if (isSkippableLine(record_text) or isFastaIdLine(record_text)) {

            // Suppress IO, we already have an ID line.
            parser_token = ParserToken::FIND_FASTA_ID_NO_READ;
            // Create a fasta record and clear the data so that a following EOF does not duplicate it.
            fasta_sequences.push_back(createFastaSequence(fasta_id_comment.first, fasta_id_comment.second, fasta_lines));
            fasta_lines.clear();
            fasta_id_comment = {};

          } else {

            // Move the fasta data line into a vector of fasta lines.
            fasta_lines.emplace_back(std::move(record_text));

          }

        }
          break;

        case ParserToken::FASTA_PARSE_ERROR: {

          // Look for the next fasta id block.
          if (isFastaIdLine(record_text)) {

            parser_token = ParserToken::FIND_FASTA_ID;

          }

        }
          break;

        case ParserToken::READ_IO_EOF: {

          // Create a fasta record if an unterminated trailing record is available.
          if (not fasta_id_comment.first.empty() and not fasta_lines.empty()) {

            fasta_sequences.push_back(createFastaSequence(fasta_id_comment.first, fasta_id_comment.second, fasta_lines));

          }

        }
          break;

      } // switch

    } while (parser_token != ParserToken::READ_IO_EOF);

  } catch(...) {

    ExecEnv::log().error("ParseFasta::readFastaFile; Unexpected IO exception");
    return false;

  }

  return true;

}


kgl::ReadFastaSequence kgl::ParseFasta::createFastaSequence( const std::string& fasta_id,
                                                             const std::string& fasta_comment,
                                                             const std::vector<std::string>& fasta_lines) {

  auto fasta_data = std::make_unique<std::string>();
  size_t fasta_size{0};

  for (const auto& line : fasta_lines) {

    fasta_size += line.size();

  }

  if (fasta_size == 0) {

    ExecEnv::log().error("ParseFasta::createFastaSequence; zero-sized fasta data");

  }

  fasta_data->reserve(fasta_size);

  for (const auto& line : fasta_lines) {

    fasta_data->append(line);

  }

  return {std::string(fasta_id), std::string(fasta_comment), std::move(fasta_data)};

}
