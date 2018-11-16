//
// Created by kellerberrin on 10/11/18.
//

#include "kgl_epigenetic_motif.h"

#include <fstream>
#include <sstream>


namespace kgl = kellerberrin::genome;


void kgl::PromoterMotif::displayTFFMotif(std::shared_ptr<const GenomeDatabase> genome_db_ptr, const std::string& file_name) {

  std::ofstream motif_file(file_name);

  for (auto contig : genome_db_ptr->getMap()) {

    TSSVector vector = contig.second->getTSSVector();

    size_t assigned_count = 0;
    for (auto tss_feature : vector) {

      ContigSize_t tss_size = (tss_feature->sequence().end()-tss_feature->sequence().begin());
      ContigOffset_t tss_offset = tss_feature->sequence().begin();

      // Check offset and size.
      if ((tss_offset + tss_size) > contig.second->sequence().length() or tss_size > contig.second->sequence().length()) {

        ExecEnv::log().warn("displayTFFMotif; contig offset: {} and region size: {} too large for contig: {} length: {}",
                            tss_offset, tss_size, contig.first, contig.second->sequence().length());
        continue;

      }

      // Get the reference DNA sequence
      std::shared_ptr<DNA5SequenceLinear> tss_sequence = contig.second->sequence().subSequence(tss_offset, tss_size);

      std::stringstream ss;

      if (tss_feature->superFeatures().size() > 0) {

        ++assigned_count;

        ss << "Gene ID: " << tss_feature->superFeatures().begin()->first;
        ss << " [" << tss_feature->superFeatures().begin()->second->sequence().begin();
        ss << ", " << tss_feature->superFeatures().begin()->second->sequence().end();
        ss << ")" << static_cast<char>(tss_feature->superFeatures().begin()->second->sequence().strand());
        ss << ", ID:" << tss_feature->id();
        ss << " [" << tss_feature->sequence().begin();
        ss << ", " << tss_feature->sequence().end();
        ss << ")" << static_cast<char>(tss_feature->sequence().strand());
        ss << " size:" << tss_size;
        long offset = tss_feature->superFeatures().begin()->second->sequence().begin();
        offset -= tss_feature->sequence().end();
        ss << " offset:" << offset;
        ss << " " << tss_sequence->getSequenceAsString();

        motif_file << ss.str() << '\n';

      }


    }

    ExecEnv::log().info("Contig: {} has {} TSS blocks defined, assigned: {}", contig.first, vector.size(), assigned_count);

  }


}


