/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KOL_PRECOMPUTED_MATRIX_TERM_SIMILARITY
#define KOL_PRECOMPUTED_MATRIX_TERM_SIMILARITY

#include <vector>
#include <exception>

#include "kol_SimilarityInterface.h"

namespace kellerberrin::ontology {

//! A class to calculate similarity between go terms for 2 sets using a precomuted term similarity matrix.
/*! \class SimilarityMatrix
	This class allows the term similarity calculation to be decoupled
	 from term set (gene) similarity measure that use them.
	 Term similarity is loaded from a matrix file.
*/
class SimilarityMatrix : public SimilarityInterface {

public:

  //! A constructor
  /*!
    Parses a matrix file and creates the SimilarityMatrix object
  */
  explicit SimilarityMatrix(const std::string &matrix_file) { readSimilarityMatrix(matrix_file); }

  ~SimilarityMatrix() override = default;

  //! A method to check if the file exists and fits the format
  /*!
    This method is used to test if a file can be used.
  */
  [[nodiscard]] bool isMatrixFileGood(const std::string &filename) const { return isFileGood(filename); }

  // termCount() == 0 is an error condition.
  [[nodiscard]] size_t termCount() const { return _matrix.size(); }
  //! A method for calculating term-to-term similarity for GO terms using a precomputed similarity matrix.
  /*!
    This method returns the term similarity as defined by the matrix.
  */
  [[nodiscard]] double calculateTermSimilarity(const std::string &goTermA, const std::string &goTermB) const override;


  //! This method projects a set of terms into it the kernel space
  /*!
    This method treats the term similarity matrix as a kernel
     and projects a set of terms into it.
  */
  [[nodiscard]] std::vector<double> projectTermSet(const std::vector<std::string> &terms) const;

private:


  OntologyMapType<std::string, std::size_t> _termToIndex;
  std::vector<std::vector<double> > _matrix;

  [[nodiscard]] bool hasTerm(const std::string &term) const {

    return _termToIndex.find(term) != _termToIndex.end();

  }

  void readSimilarityMatrix(const std::string &matrix_file);

  [[nodiscard]] bool isFileGood(const std::string &filename) const;

};

} // namespace

#endif


