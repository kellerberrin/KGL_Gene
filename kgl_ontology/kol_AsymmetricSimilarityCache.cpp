//
// Created by kellerberrin on 17/4/21.
//


#include "kel_thread_pool.h"
#include "kol_AsymmetricSimilarityCache.h"


namespace kol = kellerberrin::ontology;

//! A method for calculating term-to-term similarity for GO terms using a precomputed similarity matrix.
/*!
  This method returns the term similarity as defined by the matrix.
*/
double kol::AsymmetricSimilarityCache::calculateTermSimilarity(const std::string &row_term, const std::string &column_term) const {


  if (row_term.empty() or column_term.empty()) {

    ExecEnv::log().warn("AsymmetricSimilarityCache::calculateTermSimilarity; Empty Go term requested, Row GO term: {}, Column GO term: {}", row_term, column_term);
    return 0.0;

  }

  auto const lookup_row = row_term_to_index_.find(row_term);
  if (lookup_row == row_term_to_index_.end()) {

    ExecEnv::log().error("AsymmetricSimilarityCache::calculateTermSimilarity; Row GO term: {} not found in cached Rows.", row_term);
    if (column_term_to_index_.find(row_term) != column_term_to_index_.end()) {

      ExecEnv::log().warn("AsymmetricSimilarityCache::calculateTermSimilarity; Row GO term: {} found in Column cache.", row_term);

    }
    return 0.0;

  }

  auto const lookup_column = column_term_to_index_.find(column_term);
  if (lookup_column == column_term_to_index_.end()) {

    ExecEnv::log().error("AsymmetricSimilarityCache::calculateTermSimilarity; Column GO term: {} not found in cached Columns.", column_term);
    if (row_term_to_index_.find(column_term) != row_term_to_index_.end()) {

      ExecEnv::log().warn("AsymmetricSimilarityCache::calculateTermSimilarity; Column GO term: {} found in Row cache.", column_term);

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

double kol::AsymmetricSimilarityCache::lookUpMatrix(size_t row, size_t column) const {

  if (row >= cache_matrix_.size()) {

    ExecEnv::log().error("AsymmetricSimilarityCache::lookUpMatrix; Row Index: {} exceeds Row Size: {}", row, cache_matrix_.size());
    return 0.0;

  }

  if (column >= cache_matrix_[row]->size()) {

    ExecEnv::log().error("AsymmetricSimilarityCache::lookUpMatrix; Column Index: {} exceeds Column Size: {}", column, cache_matrix_[row]->size());
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
bool kol::AsymmetricSimilarityCache::termSimilarityCache( const std::vector<std::string>& row_terms,
                                                          const std::vector<std::string>& column_terms,
                                                          const std::shared_ptr<const TermSimilarityInterface> &term_similarity_ptr,
                                                          GO::Ontology ontology) {

  // Return the empty cache.
  if (ontology == GO::Ontology::ONTO_ERROR) {

    return false;

  }

  if (row_terms.empty() or column_terms.empty()) {

    return false;

  }

  size_t index{0};
  for (auto const &term : row_terms) {

    auto const[iter, result] = row_term_to_index_.try_emplace(term, index);
    // Check that all terms are distinct
    if (not result) {
      // Return the empty cache.
      row_term_to_index_.clear();
      return false;

    }
    ++index;

  }

  index = 0;
  for (auto const &term : column_terms) {

    auto const[iter, result] = column_term_to_index_.try_emplace(term, index);
    // Check that all terms are distinct
    if (not result) {
      // Return the empty cache.
      column_term_to_index_.clear();
      return false;

    }
    ++index;

  }




  // Create a cache matrix
  const size_t row_count = row_terms.size();
  cache_matrix_.reserve(row_count);
  // Create a pointer for the column term vector
  std::shared_ptr<const std::vector<std::string>> column_terms_ptr(std::make_shared<const std::vector<std::string>>(column_terms));
  // Default to HW threads less 1.
  ThreadPool thread_pool(ThreadPool::hardwareThreads());
  // Simplify the future vector type definition.
  using ColumnResult = std::unique_ptr<std::vector<double>>;
  using FutureResult = std::future<ColumnResult>;
  std::vector<FutureResult> future_vector;

  for (size_t row = 0; row < row_count; ++row) {

    const std::string &row_term = row_terms[row];

    FutureResult future = thread_pool.enqueueTask(&AsymmetricSimilarityCache::calcColumn,
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

std::unique_ptr<std::vector<double>> kol::AsymmetricSimilarityCache::calcColumn( const std::string& row_term,
                                                                                 const std::shared_ptr<const TermSimilarityInterface> &term_similarity_ptr,
                                                                                 const std::shared_ptr<const std::vector<std::string>>& column_terms_ptr) {

  size_t column_size = column_terms_ptr->size();
  std::unique_ptr<std::vector<double>> column_vector_ptr(std::make_unique<std::vector<double>>(column_size, 0.0));

  for (std::size_t column = 0; column < column_size; ++column) {

    const std::string &column_term = column_terms_ptr->at(column);
    double value = term_similarity_ptr->calculateNormalizedTermSimilarity(row_term, column_term);
    column_vector_ptr->at(column) = value;

  }

  return column_vector_ptr;

}
