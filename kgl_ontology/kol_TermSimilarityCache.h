//
// Created by kellerberrin on 7/4/21.
//

#ifndef KGL_KOL_TERMSIMILARITYCACHE_H
#define KGL_KOL_TERMSIMILARITYCACHE_H


#include <vector>
#include <string>
#include <iostream>
#include <condition_variable>
#include <functional>
#include <iostream>
#include <future>
#include <thread>
#include <queue>
#include <algorithm>

#include "kol_GoEnums.h"
#include "kol_GoGraph.h"
#include "kol_AnnotationData.h"
#include "kol_TermSimilarityInterface.h"

namespace kellerberrin::ontology {


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A general purpose thread safe queue for multiple consumer and producer threads.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



template<typename T> class MtQueue {

public:

  MtQueue() = default;
  ~MtQueue() = default;
  MtQueue(const MtQueue &) = delete;
  MtQueue(MtQueue &&) = delete;

  MtQueue &operator=(const MtQueue &) = delete;

  void push(T value) {

    // Scope for automatic locking/unlocking.
    {

      std::scoped_lock lock(mutex_);

      data_queue_.push(std::move(value));

      ++size_;
      ++activity_;

    }

    // The notification is sent after the queue is unlocked so that any threads on waitAndPop() can immediately execute.
    data_cond_.notify_one();


  }

  [[nodiscard]] T waitAndPop() {

    std::unique_lock<std::mutex> lock(mutex_);

    // wait on non-empty
    data_cond_.wait(lock, [this] { return not empty(); });

    T value = std::move(data_queue_.front());

    data_queue_.pop();

    --size_;
    ++activity_;

    // Unlock the mutex.
    lock.unlock();

    // Notify any other waiting threads after the queue is unlocked.
    data_cond_.notify_one();

    return value;

  }

  [[nodiscard]] bool empty() const { return size_ == 0; }

  [[nodiscard]] size_t size() const { return size_; }

  [[nodiscard]] size_t activity() const { return activity_; }

private:

  std::mutex mutex_;
  std::queue<T> data_queue_;
  std::condition_variable data_cond_;
  std::atomic<size_t> size_{0};
  std::atomic<size_t> activity_{0};

};


///////////////////////////////////////////////////////////////////////////////////////////
//
// A general purpose thread pool class.
//
///////////////////////////////////////////////////////////////////////////////////////////

class ThreadPool
{

  using Proc = std::function<void(void)>;

public:

  explicit ThreadPool(size_t threads) { startThreads(threads); }
  ~ThreadPool() noexcept { joinThreads(); }

  // Convenience routine, available hardware threads, one less than actually available.
  [[nodiscard]] static size_t hardwareThreads() { return std::thread::hardware_concurrency() - 1; }


  // Returns a std::future holding the function return value.
  template<typename F, typename... Args>
  [[nodiscard]] auto enqueueTask(F&& f, Args&&... args) -> std::future<typename std::result_of<F(Args...)>::type>
  {

    using return_type = typename std::result_of<F(Args...)>::type;
    auto task = std::make_shared<std::packaged_task<return_type()>>(std::bind(std::forward<F>(f), std::forward<Args>(args)...));
    std::future<return_type> future = task->get_future();
    work_queue_.push([task](){ (*task)(); });
    return future;

  }

  [[nodiscard]] size_t threadCount() const { return threads_.size(); }

private:

  std::vector<std::thread> threads_;
  MtQueue<Proc> work_queue_;

  void startThreads(size_t threads)
  {

    if (threads < 1) {

      threads = 1;

    }

    for(size_t i = 0; i < threads; ++i)

      threads_.emplace_back(std::thread([this]() {

        while(true)
        {

          auto workItem = work_queue_.waitAndPop();

          if (workItem == nullptr) {

            work_queue_.push(nullptr);
            break;

          }

          workItem();

        }
      }));
  }

  void joinThreads() {

    work_queue_.push(nullptr);

    for(auto& thread : threads_) {

      thread.join();

    }

  }

};


//! A multi-threaded class write a term similarity matrix to a memory cache.
/*! \class TermSimilarityCache
	This class creates a memory cache similarity matrix by calculating all similarity values
	  between pairs of terms. Warning this object can use several gigabytes of memory.
	  Matrix creation time will depend on the number of execution threads committed to matrix creation.
	  Defaults to (HW threads available - 1).
*/
class TermSimilarityCache : public TermSimilarityInterface {


public:
  //! Parameterized constructor
  /*!
    A simple parameterized constructor.
    This class takes an instance of a similarity interface.
  */
  TermSimilarityCache(const std::shared_ptr<const GoGraph> &go_graph_ptr,
                      const std::shared_ptr<const AnnotationData> &annotation_ptr,
                      const std::shared_ptr<const TermSimilarityInterface> &term_sim_ptr,
                      GO::Ontology ontology = GO::Ontology::BIOLOGICAL_PROCESS) {

    termSimilarityCache(go_graph_ptr, annotation_ptr, term_sim_ptr, ontology);

  }

  TermSimilarityCache(const TermSimilarityCache &) = delete; // Too expensive.
  ~TermSimilarityCache() override = default;

  TermSimilarityCache &operator=(const TermSimilarityCache &) = delete; // Too expensive.

  //! A method for calculating term-to-term similarity for GO terms using a precomputed similarity matrix.
  /*!
    This method returns the term similarity as defined by the matrix.
  */
  [[nodiscard]] double calculateTermSimilarity(const std::string &term_A, const std::string &term_B) const override {

    auto const lookup_A = term_to_index_.find(term_A);
    auto const lookup_B = term_to_index_.find(term_B);

    if (lookup_A != term_to_index_.end() and lookup_B != term_to_index_.end()) {

      auto const&[term_a, row] = *(lookup_A);
      auto const&[term_b, column] = *(lookup_B);

      return lookUpMatrix(row, column);

    } else {

      return 0.0;

    }

  }

  //! A method for calculating term-to-term similarity for GO terms using a precomputed similarity matrix.
  /*!
    This method returns the similarity scaled between 0 and 1 [0,1] inclusive
  */
  [[nodiscard]] double calculateNormalizedTermSimilarity(const std::string &term_A, const std::string &term_B) const override {

    return calculateTermSimilarity(term_A, term_B);

  }

  // termCount() == 0 is an error condition.
  [[nodiscard]] size_t termCount() const { return cache_matrix_.size(); }

private:

  OntologyMapType<std::string, std::size_t> term_to_index_;
  std::vector<std::unique_ptr<std::vector<double>>> cache_matrix_;

  [[nodiscard]] double lookUpMatrix(size_t row, size_t column) const {

    if (row >= cache_matrix_.size()) {

      return 0.0;

    }

    if (column >= cache_matrix_[row]->size()) {

      return 0.0;

    }

    if (row >= column) {

      std::swap(row, column);

    }

    size_t column_idx = column - row;

    return cache_matrix_[row]->at(column_idx);


  }

  //! A method to write a term similarity cache
  /*!
    Calculates square matrix of similarity terms and caches this in memory.
    Calculates the similarity between all pairs of terms.
    The complexity is O(N^2) * O(Term Pair Calculation Cost).
    O(Term Pair Calculation Cost) is usually near constant time,
     but some methods will be extremely slow and infeasible.
  */
  bool termSimilarityCache(const std::shared_ptr<const GoGraph> &go_graph_ptr,
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

    for (size_t row = 0; row < term_count; ++row) {

      const std::string &term_A = ontology_terms_ptr->at(row);

      FutureResult future = thread_pool.enqueueTask(&TermSimilarityCache::calcColumn,
                                                    row,
                                                    term_count,
                                                    term_A,
                                                    term_similarity_ptr,
                                                    ontology_terms_ptr);

      future_vector.push_back(std::move(future));

    }

    for (auto &future : future_vector) {

      auto column_ptr = future.get();
      cache_matrix_.push_back(std::move(column_ptr));

    }

    return true;

  }

  static std::unique_ptr<std::vector<double>> calcColumn(size_t row,
                                                         size_t term_count,
                                                         const std::string &term_A,
                                                         const std::shared_ptr<const TermSimilarityInterface> &term_similarity_ptr,
                                                         const std::shared_ptr<const std::vector<std::string>> &ontology_terms_ptr) {

    size_t column_size = term_count - row;
    std::unique_ptr<std::vector<double>> column_vector_ptr(std::make_unique<std::vector<double>>(column_size, 0.0));

    for (std::size_t column = row; column < term_count; ++column) {

      const std::string &term_B = ontology_terms_ptr->at(column);
      double value = term_similarity_ptr->calculateNormalizedTermSimilarity(term_A, term_B);
      size_t column_idx = column - row;
      column_vector_ptr->at(column_idx) = value;

    }

    return column_vector_ptr;

  }

};


} // namespace


#endif //TERMSIMILARITYCACHE_H
