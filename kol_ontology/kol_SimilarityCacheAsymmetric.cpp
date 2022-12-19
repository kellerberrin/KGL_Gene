//
// Created by kellerberrin on 17/4/21.
//


#include "kel_exec_env.h"

#include "../kel_thread/kel_thread_pool.h"
#include "kol_SimilarityCacheAsymmetric.h"


namespace kol = kellerberrin::ontology;

//! A method for calculating term-to-term similarity for GO terms using a precomputed similarity matrix.
/*!
  This method returns the term similarity as defined by the matrix.
*/
double kol::SimilarityCacheAsymmetric::calculateTermSimilarity(const std::string &row_term, const std::string &column_term) const {


  if (row_term.empty() or column_term.empty()) {

    ExecEnv::log().warn("SimilarityCacheAsymmetric::calculateTermSimilarity; Empty Go term requested, Row GO term: {}, Column GO term: {}", row_term, column_term);
    return 0.0;

  }

  auto const lookup_row = row_term_to_index_.find(row_term);
  if (lookup_row == row_term_to_index_.end()) {

    ExecEnv::log().error("SimilarityCacheAsymmetric::calculateTermSimilarity; Row GO term: {} not found in cached Rows.", row_term);
    if (column_term_to_index_.find(row_term) != column_term_to_index_.end()) {

      ExecEnv::log().warn("SimilarityCacheAsymmetric::calculateTermSimilarity; Row GO term: {} found in Column cache.", row_term);

    }
    return 0.0;

  }

  auto const lookup_column = column_term_to_index_.find(column_term);
  if (lookup_column == column_term_to_index_.end()) {

    ExecEnv::log().error("SimilarityCacheAsymmetric::calculateTermSimilarity; Column GO term: {} not found in cached Columns.", column_term);
    if (row_term_to_index_.find(column_term) != row_term_to_index_.end()) {

      ExecEnv::log().warn("SimilarityCacheAsymmetric::calculateTermSimilarity; Column GO term: {} found in Row cache.", column_term);

    }
    return 0.0;

  }

  auto const&[term_row, row] = *(lookup_row);
  auto const&[term_column, column] = *(lookup_column );

  return lookUpMatrix(row, column);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double kol::SimilarityCacheAsymmetric::lookUpMatrix(size_t row, size_t column) const {

  if (row >= cache_matrix_.size()) {

    ExecEnv::log().error("SimilarityCacheAsymmetric::lookUpMatrix; Row Index: {} exceeds Row Size: {}", row, cache_matrix_.size());
    return 0.0;

  }

  if (column >= cache_matrix_[row]->size()) {

    ExecEnv::log().error("SimilarityCacheAsymmetric::lookUpMatrix; Column Index: {} exceeds Column Size: {}", column, cache_matrix_[row]->size());
    return 0.0;

  }

  return cache_matrix_[row]->at(column);

}

//! A method to write a term similarity cache
/*!
  Calculates square matrix of similarity terms and caches this in memory.
  Calculates the similarity between all pairs of terms.
  This calculation is multi-threaded for performance.
  The complexity is O(Rows * Columns) * O(Term Pair Calculation Cost).
  O(Term Pair Calculation Cost) is usually near constant time,
   but some methods will be extremely slow and infeasible.
*/
bool kol::SimilarityCacheAsymmetric::termSimilarityCache(const std::vector<std::string>& row_terms,
                                                         const std::vector<std::string>& column_terms,
                                                         const std::shared_ptr<const SimilarityInterface> &term_similarity_ptr) {


  if (row_terms.empty() or column_terms.empty()) {

    ExecEnv::log().warn("SimilarityCacheAsymmetric::termSimilarityCache; empty term vectors supplied");
    return false;

  }

  size_t index{0};
  std::vector<std::string> unique_row_terms;
  for (auto const &term : row_terms) {

    auto const[iter, result] = row_term_to_index_.try_emplace(term, index);
    // Check that all terms are distinct
    if (not result) {
      // Warn and continue.
      ExecEnv::log().warn("SimilarityCacheAsymmetric::termSimilarityCache; duplicate GO term: {} in row terms vector", term);

    } else {

      unique_row_terms.push_back(term);
      ++index;

    }

  }

  index = 0;
  std::vector<std::string> unique_column_terms;
  for (auto const &term : column_terms) {

    auto const[iter, result] = column_term_to_index_.try_emplace(term, index);
    // Check that all terms are distinct
    if (not result) {
      // Warn and continue.
      ExecEnv::log().warn("SimilarityCacheAsymmetric::termSimilarityCache; duplicate GO term: {} in column terms vector", term);

    } else {

      ++index;
      unique_column_terms.push_back(term);

    }

  }




  // Create a cache matrix
  const size_t row_count = unique_row_terms.size();
  cache_matrix_.reserve(row_count);
  // Create a pointer for the column term vector
  std::shared_ptr<const std::vector<std::string>> column_terms_ptr(std::make_shared<const std::vector<std::string>>(unique_column_terms));
  // Default to HW threads less 1.
  WorkflowThreads thread_pool(WorkflowThreads::defaultThreads());
  // Simplify the future vector type definition.
  using ColumnResult = std::unique_ptr<std::vector<double>>;
  using FutureResult = std::future<ColumnResult>;
  std::vector<FutureResult> future_vector;

  for (size_t row = 0; row < row_count; ++row) {

    const std::string &row_term = unique_row_terms[row];

    FutureResult future = thread_pool.enqueueTask(&SimilarityCacheAsymmetric::calcColumn,
                                                  row_term,
                                                  term_similarity_ptr,
                                                  column_terms_ptr);

    future_vector.push_back(std::move(future));

  }

  for (auto &future : future_vector) {

    auto column_ptr = future.get();
    cache_matrix_.push_back(std::move(column_ptr));

  }

  return true;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::unique_ptr<std::vector<double>> kol::SimilarityCacheAsymmetric::calcColumn(const std::string& row_term,
                                                                                const std::shared_ptr<const SimilarityInterface> &term_similarity_ptr,
                                                                                const std::shared_ptr<const std::vector<std::string>>& column_terms_ptr) {

  size_t column_size = column_terms_ptr->size();
  std::unique_ptr<std::vector<double>> column_vector_ptr(std::make_unique<std::vector<double>>(column_size, 0.0));

  for (std::size_t column = 0; column < column_size; ++column) {

    const std::string &column_term = column_terms_ptr->at(column);
    double value = term_similarity_ptr->calculateTermSimilarity(row_term, column_term);
    column_vector_ptr->at(column) = value;

  }

  return column_vector_ptr;

}
