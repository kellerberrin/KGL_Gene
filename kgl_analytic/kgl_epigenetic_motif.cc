//
// Created by kellerberrin on 10/11/18.
//

#include "kgl_epigenetic_motif.h"

#include <fstream>
#include <sstream>


namespace kgl = kellerberrin::genome;


void kgl::PromoterMotif::displayTFFMotif(std::shared_ptr<const GenomeDatabase> genome_db_ptr, const std::string& file_name, const char delimiter) {

  std::ofstream motif_file(file_name);

  headerTFFMotif(motif_file, delimiter);

  for (auto contig : genome_db_ptr->getMap()) {

    TSSVector vector = contig.second->getTSSVector();

    size_t assigned_count = 0;
    for (auto tss_feature : vector) {

      ContigSize_t tss_size = (tss_feature->sequence().end()-tss_feature->sequence().begin());
      ContigOffset_t tss_offset = tss_feature->sequence().begin();

      // Get the reference DNA sequence
      std::shared_ptr<DNA5SequenceLinear> tss_sequence = contig.second->sequence().subSequence(tss_offset, tss_size);
      // Strand it.
      std::shared_ptr<DNA5SequenceCoding> tss_coding = tss_sequence->codingSequence(tss_feature->sequence().strand());

      if (not tss_feature->superFeatures().empty()) {

        ++assigned_count;

        motif_file << tss_feature->superFeatures().begin()->first << delimiter;
        motif_file << tss_feature->superFeatures().begin()->second->sequence().begin() << delimiter;
        motif_file << tss_feature->superFeatures().begin()->second->sequence().end() << delimiter;
        motif_file << static_cast<char>(tss_feature->superFeatures().begin()->second->sequence().strand()) << delimiter;
        long offset = tss_feature->superFeatures().begin()->second->sequence().begin();
        offset -= tss_feature->sequence().end();
        motif_file << offset << delimiter;

      } else {

        motif_file << delimiter;
        motif_file << delimiter;
        motif_file << delimiter;
        motif_file << delimiter;
        motif_file << delimiter;

      }

      motif_file << tss_feature->id() << delimiter;
      motif_file << tss_feature->sequence().begin() << delimiter;
      motif_file << tss_feature->sequence().end() << delimiter;
      motif_file << static_cast<char>(tss_feature->sequence().strand()) << delimiter;
      motif_file << tss_size << delimiter;
      motif_file << tss_coding->getSequenceAsString() << delimiter << '\n';

    }

    ExecEnv::log().info("Contig: {} has {} TSS blocks defined, assigned: {}", contig.first, vector.size(), assigned_count);

  }


}


void kgl::PromoterMotif::headerTFFMotif(std::ofstream& motif_file, const char delimiter) {

  motif_file << "GENE" << delimiter;
  motif_file << "Begin" << delimiter;
  motif_file << "End" << delimiter;
  motif_file << "Strand" << delimiter;
  motif_file << "TFFOffset" << delimiter;
  motif_file << "TFF" << delimiter;
  motif_file << "TFFBegin" << delimiter;
  motif_file << "TFFEnd" << delimiter;
  motif_file << "TFFStrand" << delimiter;
  motif_file << "TFFSize" << delimiter;
  motif_file << "TFFSequence" << delimiter << '\n';

}