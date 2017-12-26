//
// Created by kellerberrin on 26/12/17.
//

#ifndef KGL_VARIANT_FACTORY_VCF_H
#define KGL_VARIANT_FACTORY_VCF_H

#include <memory>
#include <string>
#include <map>
#include "kgl_variant_db.h"
#include "kgl_logging.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



class VcfFactory {

public:

  explicit VcfFactory();
  ~VcfFactory();

  // Functionality passed to the implmentation.
  std::shared_ptr<GenomeVariant> readParseVcf(const std::string& genome_name,
                                              std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                              const std::string& vcf_file_name,
                                              Phred_t read_quality,
                                              NucleotideReadCount_t min_read_count,
                                              double min_proportion);

private:

  class VcfFileImpl;       // Forward declaration of the VCF File reader implementation class
  std::unique_ptr<VcfFileImpl> vcf_file_impl_ptr_;    // Read VCF file PIMPL

  bool readVcfFile(const std::string& vcf_file_name);

};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_VARIANT_FACTORY_VCF_H
