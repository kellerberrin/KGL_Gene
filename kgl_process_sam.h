//
// Created by kellerberrin on 18/09/17.
//

#ifndef KGL_PROCESS_SAM_H
#define KGL_PROCESS_SAM_H

#include <memory>
#include "kgl_logging.h"
#include "kgl_read_sam.h"
#include "kgl_consume_sam.h"
#include "kgl_genome_db.h"

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


// Simple class to bolt the producer, consumer and data together.
using ContigCountData = ContigDataMap<ConsumerLocalRecord>;
using LocalConsumer = ConsumeMTSAM<ContigCountData>;
using LocalProducer = ProduceMTSAM<LocalConsumer>;

class SamCountReader {

public:

  explicit SamCountReader(Logger& log) : contig_data_ptr_(std::make_shared<ContigCountData>()) {

    consumer_ptr_ = std::shared_ptr<LocalConsumer>(std::make_shared<LocalConsumer>(log, contig_data_ptr_));
    producer_ptr_ = std::unique_ptr<LocalProducer>(std::make_unique<LocalProducer>(log, consumer_ptr_));

  }
  virtual ~SamCountReader() = default;

  std::shared_ptr<ContigCountData> readSAMFile(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                               const std::string& file_name,
                                               unsigned char readQuality = 0) {
    // Register with the genome database to setup the contig data blocks before reading.
    contig_data_ptr_->fileName(file_name);
    genome_db_ptr->registerContigData(contig_data_ptr_);
    consumer_ptr_->readQuality(readQuality);
    producer_ptr_->readSamFile(file_name);
    consumer_ptr_->finalize();
    return contig_data_ptr_;

  }

private:

  std::shared_ptr<ContigCountData> contig_data_ptr_; // Holds the raw read data.
  std::shared_ptr<LocalConsumer> consumer_ptr_; ;  // Thread safe SAM record consumer
  std::unique_ptr<LocalProducer> producer_ptr_;   //  Thread safe SAM record producer.

};

}   // namespace genome
}   // namespace kellerberrin

#endif //KGL_PROCESS_SAM_H
