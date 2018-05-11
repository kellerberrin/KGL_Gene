//
// Created by kellerberrin on 11/05/18.
//

#ifndef KGD_MT_RANDOM_H
#define KGD_MT_RANDOM_H


#include <thread>
#include <kgl_exec_env.h>
#include "kgl_mt_queue.h"
#include "kgd_random_generator.h"


namespace kellerberrin {    // organization level namespace
namespace deploid {          // project level namespace


namespace kgl = kellerberrin::genome;

// Spawn a thread to queue pre-generated random numbers in a tidal queue.
// The thread is blocked when high-tide is reached.
// The thread is re-scheduled when low_tide is reached.
template<class T>
class MTRandomGenerator : public RandomGenerator
{

public:

  MTRandomGenerator(size_t seed, size_t high_tide = HIGH_TIDE_, size_t low_tide = LOW_TIDE_);
  ~MTRandomGenerator();

  double sample();

private:

  static constexpr const size_t HIGH_TIDE_ = 10000;
  static constexpr const size_t LOW_TIDE_ = 1000;
  static constexpr const size_t BLOCK_SIZE_ = 1000;

  std::shared_ptr<T> random_generator_;
  size_t array_index_;
  size_t block_count_;
  kgl::BoundedMtQueue<std::shared_ptr<std::array<double, BLOCK_SIZE_>>> tidal_queue_;
  std::shared_ptr<std::array<double, BLOCK_SIZE_>> current_block_;
  std::shared_ptr<std::thread> producer_thread_;

  void enqueueRandom();

};

template<class T>
MTRandomGenerator<T>::MTRandomGenerator(size_t seed, size_t high_tide, size_t low_tide) : random_generator_(std::make_shared<T>(seed)),
                                                                                          array_index_(0),
                                                                                          block_count_(0),
                                                                                          tidal_queue_(high_tide, low_tide) {

  producer_thread_ = std::make_shared<std::thread>(&MTRandomGenerator<T>::enqueueRandom, this);

}

template<class T>
MTRandomGenerator<T>::~MTRandomGenerator() {

  producer_thread_->detach();

}


template<class T>
void MTRandomGenerator<T>::enqueueRandom() {

  while(true) {

    std::shared_ptr<std::array<double, BLOCK_SIZE_>>  array_ptr(std::make_shared<std::array<double, BLOCK_SIZE_>>());

    for (size_t idx = 0; idx < BLOCK_SIZE_; ++idx) {

      array_ptr->at(idx) = random_generator_->sample();

    }

    tidal_queue_.push(array_ptr);

  }

}

template<class T>
double MTRandomGenerator<T>::sample() {


  if (array_index_ == 0 or not current_block_) {

    tidal_queue_.waitAndPop(current_block_);

  }

  double rand = current_block_->at(array_index_);

  ++array_index_;

  if (array_index_ >= BLOCK_SIZE_) {

    array_index_ = 0;
    ++block_count_;
    kgl::ExecEnv::log().info("Cached Rand Block: {}, rand: {}", block_count_, rand);

  }

  return rand;

}



}   // organization level namespace
}   // project level namespace


#endif //KGD_MT_RANDOM_H
