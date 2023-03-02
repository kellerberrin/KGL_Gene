//
// Created by kellerberrin on 6/1/21.
//

#include "kgl_square_parser.h"


#include <fstream>
#include "kgl_pfgenome_aux.h"
#include "kel_exec_env.h"



namespace kgl = kellerberrin::genome;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::SquareTextIndexed::parseTextRows(const SquareTextRows& square_text) {

  // The first row is assumed to be the header.
  if (square_text.getRowVector().empty()) {

    ExecEnv::log().warn("SquareTextIndexed::parseTextRows; square text argument is empty");
    return false;

  }

  // Check field sizes are the same.
  size_t header_size = square_text.getRowVector().front().size();
  if (header_size == 0) {

    ExecEnv::log().error("SquareTextIndexed::parseTextRows; header size is zero");
    return false;

  }

  if (not square_text.checkRowSize(header_size)) {

    ExecEnv::log().error("SquareTextIndexed::parseTextRows; data row sizes different from header size: {}", header_size);
    return false;

  }

  // Index the header with field offsets.
  size_t field_index{0};
  for (auto const& header_field : square_text.getRowVector().front()) {

    auto [it, result] = square_file_header_.try_emplace(header_field, field_index);

    if (not result) {

      ExecEnv::log().error("SquareTextIndexed::parseTextRows; header field: {} is a duplicate", result);
      square_file_header_.clear();
      return false;

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

        ExecEnv::log().error("SquareTextIndexed::parseTextRows; data row: {} has index duplicate", data_row.front());
        square_file_header_.clear();
        square_indexed_rows_.clear();
        return false;

      }

    }

  }

  return true;

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

  FileDataIO file_io;
  std::shared_ptr<SquareTextRows> square_text_ptr(std::make_shared<SquareTextRows>());
  size_t counter = 0;

  if (not file_io.open(file_name)) {

    ExecEnv::log().critical("SquareTextParser::parseIndexedFile; I/O error; could not open file: {}", file_io.fileName());

  }

  while (true) {

    auto line_record = file_io.readLine();
    if (line_record.EOFRecord()) break;

    const std::string_view record_str = line_record.getView();

    if (not record_str.empty()) {

      if (record_str[0] == COMMENT_) {

        continue;  // Skip comment lines.

      }

    }

    std::vector<std::string> row_fields = Utility::charTokenizer(record_str, delimiter);
    square_text_ptr->getRowVector().push_back(std::move(row_fields));

    ++counter;

  }

  ExecEnv::log().info("Parsed: {} lines from square text file: {}", counter, file_io.fileName());

  return square_text_ptr;

}



std::optional<double> kgl::SquareTextRows::getFloat(const std::string& float_str) {


  double float_value{0.0};
  try {

    float_value = std::stod(float_str);

  } catch(...) {

    ExecEnv::log().error("SquareTextRows::getFloat; string: {}, has invalid double value", float_str);
    return std::nullopt;

  }

  return float_value;

}


std::optional<std::string> kgl::SquareTextRows::getString(const std::string& str) {

  return str;

}


std::optional<int64_t> kgl::SquareTextRows::getInteger(const std::string& int_str) {

  int64_t integer_value{0};
  try {

    integer_value = std::stoll(int_str);

  } catch(...) {

    ExecEnv::log().error("SquareTextRows::getInteger; integer string: {}, is an invalid signed integer value", int_str);
    return std::nullopt;

  }

  return integer_value;

}


std::optional<size_t> kgl::SquareTextRows::getSize(const std::string& size_str) {

  size_t size_value{0};
  try {

    size_value = std::stoull(size_str);

  } catch(...) {

    ExecEnv::log().error("SquareTextRows::getSize; size string: {}, is an invalid unsigned integer value", size_str);
    return std::nullopt;

  }

  return size_value;

}

