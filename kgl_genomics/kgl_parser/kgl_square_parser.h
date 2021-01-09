//
// Created by kellerberrin on 6/1/21.
//

#ifndef KGL_SQUARE_PARSER_H
#define KGL_SQUARE_PARSER_H




#include <string>
#include <vector>
#include <map>

#include "kgl_data_file_impl.h"


namespace kellerberrin::genome {   //  organization::project level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// General purpose parser for a small square text files.
// The file format is as follows:
// The first line is a delimited text header.
// Subsequent rows have the same number of fields as the header.
// The first column of the file is assumed to be a record index and is used to create an indexed map of fields.
// Example:
// Header (1st) line: "index, field1, field2, field3"
// First data line: "idx1, value1, value2, value3"
//
// Comments can be insert in the file by starting the line with the character '#', e.g. "# This is a comment\n"
// Comments are ignored by the parser.
//
// Note that the file fields will generally be comma ',' (.csv) or tab '\t' (.tsv) delimited.
// Leading and trailing (but not internal) whitespace will be stripped off fields.
// Rows are accessed by index.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using SquareRowVector = std::vector<std::vector<std::string>>;

class SquareTextRows {

public:

  SquareTextRows() = default;
  ~SquareTextRows() = default;

  [[nodiscard]] const SquareRowVector& getRowVector() const { return square_rows_; }
  [[nodiscard]] SquareRowVector& getRowVector() { return square_rows_; }
  [[nodiscard]] bool checkRowSize(size_t size) const;   // Check that all rows have a constant size.

private:

  SquareRowVector square_rows_;

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// An indexed structure with header.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using SquareFileHeader = std::map<std::string, size_t>;
using SquareIndexedRows = std::map<std::string, std::vector<std::string>>;


class SquareTextIndexed {

public:

  SquareTextIndexed(const SquareTextRows& square_text);
  ~SquareTextIndexed() = default;

  [[nodiscard]] const SquareFileHeader& getHeaderMap() const { return square_file_header_; }
  [[nodiscard]] const SquareIndexedRows& getRowMap() const { return square_indexed_rows_; }
  [[nodiscard]] std::optional<size_t> fieldOffset(const std::string& field_header) const;
  [[nodiscard]] const std::vector<std::string>& getRow(const std::string& row_index) const;


private:

  SquareFileHeader square_file_header_;  // Always assumed to be the first line.
  SquareIndexedRows square_indexed_rows_;

};



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The actual base class parser.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class SquareTextParser {

public:

  SquareTextParser() = default;
  ~SquareTextParser() = default;

  [[nodiscard]] std::shared_ptr<SquareTextRows> parseFlatFile(const std::string& file_name, char delimiter);

  constexpr static const char DELIMITER_CSV{','};
  constexpr static const char DELIMITER_TSV{'\t'};

private:

  FileDataIO file_io_;
  constexpr static const char COMMENT_{'#'}; // If first character then line ignored.


};



}   // end namespace




#endif //KGL_SQUARE_PARSER_H
