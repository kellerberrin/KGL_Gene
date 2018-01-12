//
// Created by kellerberrin on 26/12/17.
//

#ifndef KGL_VARIANT_FACTORY_BAM_H
#define KGL_VARIANT_FACTORY_BAM_H


#include <memory>
#include <string>
#include <map>
#include "kgl_variant_db.h"
#include "kgl_logging.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

class BamFactory {

public:

  explicit BamFactory() = default;
  ~BamFactory() = default;

  // Functionality passed to the implmentation.
  std::shared_ptr<GenomeVariant> readParseBam(const std::string& genome_name,
                                              std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                              const std::string& bam_file_name,
                                              Phred_t read_quality,
                                              Phred_t variant_quality,
                                              NucleotideReadCount_t min_read_count,
                                              double min_proportion);

private:


};


}   // namespace genome
}   // namespace kellerberrin


#endif // KGL_VARIANT_FACTORY_BAM_H
