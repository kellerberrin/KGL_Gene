//
// Created by kellerberrin on 28/09/22.
//

#include "kel_utility.h"
#include "kgl_fasta.h"


namespace kgl = kellerberrin::genome;


std::shared_ptr<kgl::GenomeReference> kgl::ParseFasta::readFastaFile( const std::string& organism, const std::string& fasta_file_name) {

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


bool kgl::ParseFasta::writeFastaFile( const std::string& fasta_file_name,
                                         const std::vector<WriteFastaSequence>& fasta_sequences) {

  return true;

}

bool kgl::ParseFasta::readFastaFile( const std::string& fasta_file_name,
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


kgl::ReadFastaSequence kgl::ParseFasta::createFastaSequence( const std::string& fasta_id,
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

