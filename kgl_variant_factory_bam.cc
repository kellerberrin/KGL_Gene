//
// Created by kellerberrin on 26/12/17.
//



#include "kgl_utility.h"
#include "kgl_sam_process.h"
#include "kgl_variant_factory_bam.h"
#include "kgl_variant_single.h"
#include "kgl_variant_factory_single.h"


namespace kgl = kellerberrin::genome;


// Functionality passed to the implmentation.
std::shared_ptr<kgl::GenomeVariant> kgl::BamFactory::readParseBam(const std::string& genome_name,
                                                                  std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                                                  const std::string& bam_file_name,
                                                                  Phred_t read_quality,
                                                                  NucleotideReadCount_t min_read_count,
                                                                  double min_proportion) {

  std::shared_ptr<GenomeVariant> genome_single_variants = kgl::GenomeVariant::emptyGenomeVariant(genome_name, genome_db_ptr);

  // To be implemented.
  // To be implemented.
  ExecEnv::log().warn("BAM file factory is not yet implmented, 0 variants will be returned for BAM file :{}", bam_file_name);

  ExecEnv::log().info("Bam file {} has: {} raw variants", bam_file_name, genome_single_variants->size());
  return genome_single_variants;

}
