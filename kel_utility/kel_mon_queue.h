//
// Created by kellerberrin on 18/12/20.
//

#ifndef KEL_MON_QUEUE_H
#define KEL_MON_QUEUE_H

#include <string>

#include "kel_mt_queue.h"

namespace kellerberrin {   //  organization level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Simple object collects tidal queue statistics to facilitate fine-tuning thread allocation in the parser io chains.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
class MonitorTidalQueue {

public:

  MonitorTidalQueue(const T& queue) : high_tide_(queue.highTide()),
                                      low_tide_(queue.lowTide()),
                                      queue_name_(queue.queueName()) {}
  ~MonitorTidalQueue() { if (queueSamples() > 0) displayQueueStats(); }  // Display queue stats on destruction.

  [[nodiscard]] size_t highTide() const { return high_tide_; }
  [[nodiscard]] size_t lowTide() const { return low_tide_; }
  [[nodiscard]] const std::string& queueName() const { return queue_name_; }
  [[nodiscard]] size_t cumulativeQueueSize() const { return cumulative_queue_size_; }
  [[nodiscard]] size_t queueSamples() const { return queue_samples_; }
  [[nodiscard]] double averageSize() const {

    if (queue_samples_ > 0) {

      return static_cast<double>(cumulative_queue_size_) / static_cast<double>(queue_samples_);

    }

    return 0.0;

  }

  void sampleQueue(const T& queue) {

    cumulative_queue_size_ += queue.size();
    ++queue_samples_;

  }

  void displayQueueStats() const {

    ExecEnv::log().info("Monitor queue: {}, high tide: {}, low tide: {}, samples: {}, average size: {}",
                        queueName(), highTide(), lowTide(), queueSamples(), averageSize());

  }

private:

  const size_t high_tide_;
  const size_t low_tide_;
  const std::string queue_name_;
  size_t cumulative_queue_size_{0};
  size_t queue_samples_{0};

};


} // namespace


#endif //KEL_MON_QUEUE_H
