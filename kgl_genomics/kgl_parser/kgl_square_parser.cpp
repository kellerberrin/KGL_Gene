//
// Created by kellerberrin on 6/1/21.
//

#include "kgl_square_parser.h"


#include <fstream>
#include "kgl_genome_aux_csv.h"
#include "kel_exec_env.h"



namespace kgl = kellerberrin::genome;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

kgl::SquareTextIndexed::SquareTextIndexed(const SquareTextRows& square_text) {

// The first row is assumed to be the header.
  if (square_text.getRowVector().empty()) {

    ExecEnv::log().warn("SquareTextIndexed::SquareTextIndexed; square text argument is empty");
    return;

  }

  // Check field sizes are the same.
  size_t header_size = square_text.getRowVector().front().size();
  if (header_size == 0) {

    ExecEnv::log().error("SquareTextIndexed::SquareTextIndexed; header size is zero");
    return;

  }

  if (not square_text.checkRowSize(header_size)) {

    ExecEnv::log().error("SquareTextIndexed::SquareTextIndexed; data row sizes different from header size: {}", header_size);
    return;

  }

  // Index the header with field offsets.
  size_t field_index{0};
  for (auto const& header_field : square_text.getRowVector().front()) {

    auto [it, result] = square_file_header_.try_emplace(header_field, field_index);

    if (not result) {

      ExecEnv::log().error("SquareTextIndexed::SquareTextIndexed; header field: {} is a duplicate", result);
      square_file_header_.clear();
      return;

    }

  }

  // Data rows are assumed to be indexed by the first column
  bool first_row = true;
  for (auto const& data_row : square_text.getRowVector()) {

    if (first_row) {

      first_row = false;

    } else {

      auto [it, result] = square_indexed_rows_.try_emplace(data_row.front(), data_row);

      if (not result) {

        ExecEnv::log().error("SquareTextIndexed::SquareTextIndexed; data row: {} has index duplicate", data_row.front());
        square_file_header_.clear();
        square_indexed_rows_.clear();
        return;

      }

    }

  }

}


// Returns the field offset
std::optional<size_t> kgl::SquareTextIndexed::fieldOffset(const std::string& field_name) const {

  auto result = square_file_header_.find(field_name);

  if (result != square_file_header_.end()) {

    return result->second;

  }

  ExecEnv::log().warn("SquareTextIndexed::fieldOffset; could not find field: {}", field_name);

  return std::nullopt;

}


const std::vector<std::string>& kgl::SquareTextIndexed::getRow(const std::string& row_index) const {

  auto result = square_indexed_rows_.find(row_index);

  if (result != square_indexed_rows_.end()) {

    return result->second;

  }

  ExecEnv::log().warn("SquareTextIndexed::getRow; could not find row index: {}", row_index);

  static std::vector<std::string> empty_row;

  return empty_row;


}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Check that all rows have a constant size.
[[nodiscard]] bool kgl::SquareTextRows::checkRowSize(size_t size) const {

  for (auto const& row : square_rows_) {

    if (size != row.size()) return false;

  }

  return true;

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::shared_ptr<kgl::SquareTextRows> kgl::SquareTextParser::parseFlatFile(const std::string& file_name, char delimiter) {

  std::shared_ptr<SquareTextRows> square_text_ptr(std::make_shared<SquareTextRows>());
  size_t counter = 0;

  if (not file_io_.commenceIO(file_name)) {

    ExecEnv::log().critical("SquareTextParser::parseIndexedFile; I/O error; could not open file: {}", file_io_.fileName());

  }

  while (true) {

    auto line_record = file_io_.readIORecord();
    if (not line_record) break;

    auto const& [line_count, line_ptr] = line_record.value();

    const std::string& record_str = *line_ptr;

    if (not record_str.empty()) {

      if (record_str[0] == COMMENT_) {

        continue;  // Skip comment lines.

      }

    }

    std::vector<std::string> row_fields = Utility::char_tokenizer(record_str, delimiter);
    square_text_ptr->getRowVector().push_back(std::move(row_fields));

    ++counter;

  }

  ExecEnv::log().info("Parsed: {} lines from square text file: {}", counter, file_io_.fileName());

  return square_text_ptr;

}




