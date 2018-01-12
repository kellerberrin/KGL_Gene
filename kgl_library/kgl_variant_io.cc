//
// Created by kellerberrin on 12/11/17.
//


#include <memory>
#include <fstream>
#include "kgl_patterns.h"
#include "kgl_variant_compound.h"
#include "kgl_variant_db.h"
#include "kgl_filter.h"
#include "kgl_gff_fasta.h"

namespace kgl = kellerberrin::genome;


bool kgl::GenomeVariant::outputCSV(const std::string& file_name, VariantOutputIndex output_index, bool detail) const {

  // open the file.
  std::fstream out_file(file_name, std::fstream::out | std::fstream::app);
  if (!out_file) {

    ExecEnv::log().error("Cannot open output CSV file (--outCSVFile): {}", file_name);
    return false;

  }

  out_file << output(',', output_index, detail);

  return out_file.good();

}


std::ostream& operator<<(std::ostream &os, const kgl::GenomeVariant& genome_variant) {

  os << "(Variants are displayed with the start offset = 1 convention)\n";
  os << genome_variant.output(' ', kgl::VariantOutputIndex::START_1_BASED, false);
  os.flush();

  return os;

}

std::ostream & operator<<(std::ostream &os, std::shared_ptr<const kellerberrin::genome::GenomeVariant> genome_variant_ptr) {

  return operator<<(os, *genome_variant_ptr);

}

std::ostream & operator<<(std::ostream &os, std::shared_ptr<kellerberrin::genome::GenomeVariant> genome_variant_ptr) {

  return operator<<(os, *genome_variant_ptr);

}
