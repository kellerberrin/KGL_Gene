//
// Created by kellerberrin on 19/4/21.
//


#include "kol_OntologyTypes.h"
#include "kol_SimilarityMatrix.h"

#include <fstream>
#include <boost/tokenizer.hpp>


namespace kol = kellerberrin::ontology;


//! A method for calculating term-to-term similarity for GO terms using a precomputed similarity matrix.
/*!
  This method returns the term similarity as defined by the matrix.
*/
double kol::SimilarityMatrix::calculateTermSimilarity(const std::string &goTermA, const std::string &goTermB) const {

  if (hasTerm(goTermA) && hasTerm(goTermB)) {

    auto const&[aterm, row] = *(_termToIndex.find(goTermA));
    auto const&[bterm, column] = *(_termToIndex.find(goTermB));

    return _matrix[row][column];

  } else {

    return 0.0;

  }

}


//! This method projects a set of terms into it the kernel space
/*!
  This method treats the term similarity matrix as a kernel
   and projects a set of terms into it.
*/
std::vector<double> kol::SimilarityMatrix::projectTermSet(const std::vector<std::string> &terms) const {

  std::vector<double> projectedVec(_matrix.size(), 0.0);

  for (std::size_t col = 0; col < _matrix.size(); ++col) {
    double sum = 0.0;

    for (auto const &prot : terms) {

      if (_termToIndex.find(prot) != _termToIndex.end()) {

        auto const&[term, row] = *(_termToIndex.find(prot));
        sum += _matrix.at(row).at(col);

      }

    }

    projectedVec[col] = sum;

  }
  return projectedVec;
}



void kol::SimilarityMatrix::readSimilarityMatrix(const std::string &matrix_file) {

  std::ifstream in(matrix_file.c_str());

  // Check the file stream.
  if (not in.good()) {

    return;

  }
  //Tokenizer type
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  boost::char_separator<char> tab_sep("\t", "", boost::keep_empty_tokens);
  tokenizer::iterator it;
  std::string line;

  std::size_t line_count = 0;

  //first line
  getline(in, line);
  tokenizer tokens(line, tab_sep);
  it = tokens.begin();

  std::string firstTerm = *it;
  _termToIndex[firstTerm] = line_count;
  ++line_count;
  ++it;

  std::vector<double> firstRow;
  while (it != tokens.end()) {

    std::string strval = *it;
    double val = atof(strval.c_str());
    firstRow.push_back(val);
    ++it;

  }

  std::size_t nRows = firstRow.size();
  _matrix.reserve(nRows);
  _matrix.push_back(firstRow);

  while (in.good()) {

    getline(in, line);
    tokenizer tokens(line, tab_sep);
    it = tokens.begin();
    if (it == tokens.end()) { continue; }

    std::string term = *it;
    _termToIndex[term] = line_count;
    ++line_count;
    ++it;

    std::vector<double> row;
    row.reserve(nRows);
    while (it != tokens.end()) {

      std::string strval = *it;
      double val = atof(strval.c_str());
      row.push_back(val);
      ++it;

    }

    _matrix.push_back(row);

  }
  //std::cout << "rows " << cache_matrix_.size() << std::endl;
  //std::cout << "cols " << cache_matrix_.at(0).size() << std::endl;
  //std::cout << "terms " << term_to_index_.size() << std::endl;
}


bool kol::SimilarityMatrix::isFileGood(const std::string &filename) const {

  std::ifstream in(filename.c_str());
  if (not in.good()) {

    return false;

  }

  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  boost::char_separator<char> tab_sep("\t", "", boost::keep_empty_tokens);
  tokenizer::iterator it;
  std::string line;

  std::size_t line_count = 0;

  //first line
  getline(in, line);
  tokenizer tokens(line, tab_sep);
  it = tokens.begin();

  std::string firstTerm = *it;
  ++line_count;
  ++it;

  std::vector<double> firstRow;
  while (it != tokens.end()) {
    std::string strval = *it;
    double val = atof(strval.c_str());
    firstRow.push_back(val);
    ++it;
  }

  try {

    while (in.good() && line_count < 5) {

      getline(in, line);
      tokenizer tokens(line, tab_sep);
      it = tokens.begin();
      if (it == tokens.end()) { continue; }

      std::string term = *it;
      ++line_count;
      ++it;

      std::vector<double> row;
      row.reserve(firstRow.size());

      while (it != tokens.end()) {
        std::string strval = *it;
        double val = atof(strval.c_str());
        row.push_back(val);
        ++it;
      }

      if (row.size() != firstRow.size()) {

        return false;

      }

    }

  }
  catch (std::exception &e) {

    return false;

  }

  return true;
}
