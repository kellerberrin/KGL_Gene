//
// Created by kellerberrin on 17/4/21.
//

#include "kel_exec_env.h"
#include "kel_thread_pool.h"
#include "kol_TermSimilarityCache.h"


namespace kol = kellerberrin::ontology;

//! A method for calculating term-to-term similarity for GO terms using a precomputed similarity matrix.
/*!
  This method returns the term similarity as defined by the matrix.
*/
double kol::TermSimilarityCache::calculateTermSimilarity(const std::string &row_term, const std::string &column_term) const {


  auto const lookup_row = term_to_index_.find(row_term);
  if (lookup_row == term_to_index_.end()) {

    ExecEnv::log().error("TermSimilarityCache::calculateTermSimilarity; Row GO: term not found: {}", row_term);
    return 0.0;

  }

  auto const lookup_column = term_to_index_.find(column_term);
  if (lookup_column == term_to_index_.end()) {

    ExecEnv::log().error("TermSimilarityCache::calculateTermSimilarity; Column GO: term not found: {}", column_term);
    return 0.0;

  }

  auto const&[term_row, row] = *(lookup_row);
  auto const&[term_column, column] = *(lookup_column);

  return lookUpMatrix(row, column);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double kol::TermSimilarityCache::lookUpMatrix(size_t row, size_t column) const {


  size_t adj_row;
  size_t adj_column;
  if (column >= row) {

    adj_row = column;
    adj_column = row;

  } else {

    adj_row = row;
    adj_column = column;

  }

  size_t row_idx = adj_row - adj_column;

  if (adj_column >= cache_matrix_.size()) {

    ExecEnv::log().error("TermSimilarityCache::lookUpMatrix; Column: {}, Column Index {}, exceeds Column Size: {}",
                         column, adj_column, cache_matrix_.size());
    return 0.0;

  }

  if (row_idx >= cache_matrix_[adj_column]->size()) {

    ExecEnv::log().error("TermSimilarityCache::lookUpMatrix; Row: {}, Row Index {}, exceeds Row Size: {}",
                         row, row_idx, cache_matrix_[adj_column]->size());
    return 0.0;

  }

  return cache_matrix_[adj_column]->at(row_idx);

}

//! A method to write a term similarity cache
/*!
  Calculates square matrix of similarity terms and caches this in memory.
  Calculates the similarity between all pairs of terms.
  The complexity is O(N^2) * O(Term Pair Calculation Cost).
  O(Term Pair Calculation Cost) is usually near constant time,
   but some methods will be extremely slow and infeasible.
*/
bool kol::TermSimilarityCache::termSimilarityCache(const std::shared_ptr<const GoGraph> &go_graph_ptr,
                                                   const std::shared_ptr<const AnnotationData> &annotation_ptr,
                                                   const std::shared_ptr<const TermSimilarityInterface> &term_similarity_ptr,
                                                   GO::Ontology ontology) {

  // Return the empty cache.
  if (ontology == GO::Ontology::ONTO_ERROR) {

    return false;

  }

  auto const ontology_terms_ptr(std::make_shared<const std::vector<std::string>>(annotation_ptr->getOntologyTerms(*go_graph_ptr, ontology)));

  size_t index{0};
  for (auto const &term : *ontology_terms_ptr) {

    auto const[iter, result] = term_to_index_.try_emplace(term, index);
    // Check that all terms are distinct
    if (not result) {
      // Return the empty cache.
      term_to_index_.clear();
      return false;

    }
    ++index;

  }

  // Create a cache matrix
  const size_t term_count = ontology_terms_ptr->size();
  cache_matrix_.reserve(term_count);
  // Default to HW threads less 1.
  ThreadPool thread_pool(ThreadPool::hardwareThreads());
  // Simplify the future vector type definition.
  using ColumnResult = std::unique_ptr<std::vector<double>>;
  using FutureResult = std::future<ColumnResult>;
  std::vector<FutureResult> future_vector;

  for (size_t column = 0; column < term_count; ++column) {

    const std::string &column_term = ontology_terms_ptr->at(column);

    FutureResult future = thread_pool.enqueueTask(&TermSimilarityCache::calcColumn,
                                                  column,
                                                  column_term,
                                                  term_similarity_ptr,
                                                  ontology_terms_ptr);

    future_vector.push_back(std::move(future));

  }

  for (auto &future : future_vector) {

    auto row_ptr = future.get();
    cache_matrix_.push_back(std::move(row_ptr));

  }

  return true;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::unique_ptr<std::vector<double>> kol::TermSimilarityCache::calcColumn( size_t column,
                                                                           const std::string &column_term,
                                                                           const std::shared_ptr<const TermSimilarityInterface> &term_similarity_ptr,
                                                                           const std::shared_ptr<const std::vector<std::string>> &ontology_terms_ptr) {

  size_t term_count = ontology_terms_ptr->size();
  size_t row_size = term_count - column;
  std::unique_ptr<std::vector<double>> row_vector_ptr(std::make_unique<std::vector<double>>(row_size, 0.0));

  for (std::size_t row = 0; row < row_size; ++row) {

    const std::string &row_term = ontology_terms_ptr->at(row + column);
    double value = term_similarity_ptr->calculateTermSimilarity(row_term, column_term);
    row_vector_ptr->at(row) = value;

  }

  return row_vector_ptr;

}
