//
// Created by kellerberrin on 29/6/20.
//

#ifndef KGL_VARIANT_FACTORY_POPULATION_H
#define KGL_VARIANT_FACTORY_POPULATION_H


#include "../../kel_thread/kel_queue_mt_safe.h"
#include "kgl_variant.h"
#include "kgl_variant_db_population.h"

#include <memory>
#include <string>
#include <vector>
#include <optional>

namespace kellerberrin::genome {   //  organization::project level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This object de-queues variants from the parser threads and indexes them into a population structure.
// This indexing is done using a single thread.
// This avoids the parser threads stalling on the indexing resource constraint.

using QueuedVariant = std::optional<std::unique_ptr<const Variant>>;

class IndexVariants {

public:

  explicit IndexVariants(const std::shared_ptr<UnphasedPopulation> population_ptr)
  : population_ptr_(population_ptr), variant_queue_(HIGH_TIDE_, LOW_TIDE_) {}
  ~IndexVariants() = default;

  // Spawns indexing thread.
  void commenceIndexing() ;

  // filter are enqueued by the parser threads.
  void enqueueVariant(std::unique_ptr<const Variant>&& variant_ptr);

  // After the parser threads return, then halt and join.
  void haltProcessingAndJoin();

private:

  const std::shared_ptr<UnphasedPopulation> population_ptr_;
  QueueTidal<QueuedVariant> variant_queue_; // Parsed variant queue
  std::unique_ptr<std::thread> index_thread_ptr_;


  static constexpr const long HIGH_TIDE_{1000000};          // Maximum QueueTidal size
  static constexpr const long LOW_TIDE_{100000};            // Low water mark to begin queueing VCF records
  // Placed onto the queue by the producers to signal the end of processing.
  // Must be one of these for each producer thread.
  static constexpr const std::nullopt_t NULL_TOKEN{std::nullopt};

  // Remove a variant from the queue for indexing.
  [[nodiscard]] QueuedVariant readVariant() { return variant_queue_.waitAndPop(); }
  // Index variant code is single threaded.
  void indexVariants();

};



} // end namespace


#endif //KGL_KGL_VARIANT_FACTORY_POPULATION_H
