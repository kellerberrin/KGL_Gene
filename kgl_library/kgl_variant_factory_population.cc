//
// Created by kellerberrin on 29/6/20.
//

#include "kgl_variant_factory_population.h"


namespace kgl = kellerberrin::genome;

void kgl::IndexVariants::commenceIndexing() {

  // Commence indexing with just one thread.
  index_thread_ptr_ = std::make_unique<std::thread>(&IndexVariants::indexVariants, this);


}


void kgl::IndexVariants::haltProcessingAndJoin() {

  // Shut down the thread.
  variant_queue_.push(NULL_TOKEN);

  // Join and proceed.
  index_thread_ptr_->join();

}


void kgl::IndexVariants::enqueueVariant(std::unique_ptr<const Variant>&& variant_ptr) {

  variant_queue_.push(std::move(variant_ptr));

}



void kgl::IndexVariants::indexVariants() {

  while (true) {

    QueuedVariant variant_opt = readVariant();

    // check for NULL condition.
    if (not variant_opt) {

        break;

    }

    // Index the variant.
 //   if (not population_ptr_->addVariant(std::move(variant_opt.value()))) {

//      ExecEnv::log().error("IndexVariants::indexVariants, problem adding variant to population.");

//    }

  } // end while

}

