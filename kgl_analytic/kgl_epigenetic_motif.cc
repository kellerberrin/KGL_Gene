//
// Created by kellerberrin on 10/11/18.
//

#include "kgl_epigenetic_motif.h"


namespace kgl = kellerberrin::genome;


void kgl::PromoterMotif::displayTFFMotif(std::shared_ptr<const GenomeDatabase> genome_db_ptr) {


  for (auto contig : genome_db_ptr->getMap()) {

    TSSVector vector = contig.second->getTSSVector();

    size_t assigned_count = 0;
    for (auto tss_feature : vector) {

      if (tss_feature->superFeatures().size() > 0) {

        ++assigned_count;
        long offset = tss_feature->superFeatures().begin()->second->sequence().begin();
        offset -= tss_feature->sequence().end();
        ExecEnv::log().info("Gene ID: {} [{},{}){}, ID: {} [{},{}){} size: {} offset: {}",
                            tss_feature->superFeatures().begin()->first,
                            tss_feature->superFeatures().begin()->second->sequence().begin(),
                            tss_feature->superFeatures().begin()->second->sequence().end(),
                            static_cast<char>(tss_feature->superFeatures().begin()->second->sequence().strand()),
                            tss_feature->id(),
                            tss_feature->sequence().begin(),
                            tss_feature->sequence().end(),
                            static_cast<char>(tss_feature->sequence().strand()),
                            (tss_feature->sequence().end()-tss_feature->sequence().begin()),
                            offset);

      }

    }

    ExecEnv::log().info("Contig: {} has {} TSS blocks defined, assigned: {}", contig.first, vector.size(), assigned_count);

  }


}


