//
// Created by kellerberrin on 19/4/21.
//

#include "kol_SimilarityWriter.h"

namespace kol = kellerberrin::ontology;

//! A method to write a term similarity matrix
/*!
  Calculates and writes the similarity matrix to file.
  Calculates the similarity between all pairs of terms.
  The complexity is O(N^2) * O(Term Pair Calculation Cost).
  O(Term Pair Calculation Cost) is usually near constant time,
   but some method will be extremely slow and require other methods
*/


void kol::SimilarityWriter::writeSimilarityMatrix(const std::shared_ptr<const SimilarityInterface> &termSim,
                                                  const std::string &fileName,
                                                  GO::Ontology ontology) const {

  std::vector<std::string> ontologyTerms = annotation_ptr_->getOntologyTerms(*graph_ptr_, ontology);

  // Initialze a matrix
  std::size_t nTerms = ontologyTerms.size();
  std::vector<std::vector<double> > matrix;
  matrix.reserve(nTerms);

  for (std::size_t i = 0; i < nTerms; ++i) {

    matrix.emplace_back(std::vector<double>());

    for (std::size_t j = 0; j < nTerms; ++j) {

      matrix[i].push_back(0.0);

    }

  }

  // --- Main Calculation ---
  for (std::size_t i = 0; i < nTerms; ++i) {

    std::string termA = ontologyTerms.at(i);
    matrix[i][i] = termSim->calculateTermSimilarity(termA, termA);

    for (std::size_t j = i + 1; j < nTerms; ++j) {

      std::string termB = ontologyTerms.at(j);
      double value = termSim->calculateTermSimilarity(termA, termB);
      //cout << termA << "\t" << termB << "\t" << value << endl;
      matrix[i][j] = value;
      matrix[j][i] = value;

    }

  }

  // write the matrix
  writeMatrix(fileName, matrix, ontologyTerms);

}

//! A method to write a matrix to file
/*!
  writes the matrix to file.
*/

void kol::SimilarityWriter::writeMatrix(const std::string &fname,
                                        const std::vector<std::vector<double> > &matrix,
                                        const std::vector<std::string> &terms) {

  std::ofstream out(fname.c_str());

  std::size_t matrix_size = matrix.size();
  for (std::size_t i = 0; i < matrix_size; ++i) {

    std::string termA = terms.at(i);
    out << termA;
    for (std::size_t j = 0; j < matrix_size; ++j) {

      out << '\t' << matrix[i][j];

    }

    out << '\n';

  }

}
