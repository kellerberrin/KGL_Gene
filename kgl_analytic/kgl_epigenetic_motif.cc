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
    for (auto const tss_feature : vector) {

      ContigSize_t tss_size = (tss_feature->sequence().end()-tss_feature->sequence().begin());
      ContigOffset_t tss_offset = tss_feature->sequence().begin();

      // Get the motif DNA sequence.
      std::shared_ptr<const DNA5SequenceLinear> tss_sequence = contig.second->sequence().subSequence(tss_offset, tss_size);
      // Strand it.
      std::shared_ptr<const DNA5SequenceCoding> tss_coding = tss_sequence->codingSequence(tss_feature->sequence().strand());

#define PREFACE_SIZE 10

      ContigOffset_t preface_offset;
      if (tss_feature->sequence().strand() == StrandSense::REVERSE) {

        preface_offset = tss_feature->sequence().end();

      } else {

        preface_offset = tss_feature->sequence().begin() -  PREFACE_SIZE;

      }

      // Get the preface DNA sequence
      std::shared_ptr<const DNA5SequenceLinear> tss_preface = contig.second->sequence().subSequence(preface_offset, PREFACE_SIZE);
      // Strand it.
      std::shared_ptr<const DNA5SequenceCoding> tss_preface_coding = tss_preface->codingSequence(tss_feature->sequence().strand());


      if (not tss_feature->superFeatures().empty()) {

        ++assigned_count;

        std::shared_ptr<const Feature> gene = tss_feature->superFeatures().begin()->second;

        if (not gene->isGene()) {

          ExecEnv::log().warn("PromoterMotif::displayTFFMotif; Superfeature: {} Is Not a GENE for TSS motif: {}", gene->id(), tss_feature->id());

        }

        std::shared_ptr<const OntologyRecord> ontology_record;
        std::string symbolicReference;
        std::string description;
        std::string altSymbolicReference;

        if (not genome_db_ptr->geneOntology().getGafFeatureVector(gene->id(), ontology_record)) {

          ExecEnv::log().vwarn("PromoterMotif::displayTFFMotif; No Ontology record for GENE {}", gene->id());

        } else if (ontology_record) {

          // Strip out delimiter chars
          std::string search = { delimiter };
          symbolicReference = ontology_record->symbolicReference();
          symbolicReference = Utility::findAndReplaceAll(symbolicReference, search, ";");
          description = ontology_record->description();
          description = Utility::findAndReplaceAll(description, search, ";");
          altSymbolicReference = ontology_record->altSymbolicReference();
          altSymbolicReference = Utility::findAndReplaceAll(altSymbolicReference, search, ";");

        }

        motif_file << gene->id() << delimiter;
        motif_file << description << delimiter;
        motif_file << symbolicReference << delimiter;
        motif_file << altSymbolicReference << delimiter;
        motif_file << gene->sequence().begin() << delimiter;
        motif_file << gene->sequence().end() << delimiter;
        motif_file << static_cast<char>(gene->sequence().strand()) << delimiter;

        if (gene->sequence().strand() != tss_feature->sequence().strand()) {

          ExecEnv::log().warn("PromoterMotif::displayTFFMotif; Superfeature (Gene): {} has different strand to TSS motif: {}",
                              gene->id(), tss_feature->id());

        }

        long prime_5_offset;
        if (tss_feature->sequence().strand() == StrandSense::REVERSE) {

          prime_5_offset = tss_feature->sequence().end();
          prime_5_offset -= gene->sequence().end();

        } else {

          prime_5_offset = gene->sequence().begin();
          prime_5_offset -= tss_feature->sequence().begin();

        }

        motif_file << prime_5_offset << delimiter;

        long offset;
        offset = gene->sequence().begin();
        offset -= tss_feature->sequence().begin();

        motif_file << offset << delimiter;


      } else {

        motif_file << delimiter; // gene id
        motif_file << delimiter; // description
        motif_file << delimiter; // symbolic
        motif_file << delimiter; // altsymbolic
        motif_file << delimiter; // gene begin
        motif_file << delimiter; // gene end
        motif_file << delimiter; // gene strand
        motif_file << delimiter; // 5 prime distance
        motif_file << delimiter; // Offset

      }

      motif_file << tss_feature->id() << delimiter;
      motif_file << tss_feature->sequence().begin() << delimiter;
      motif_file << tss_feature->sequence().end() << delimiter;
      motif_file << static_cast<char>(tss_feature->sequence().strand()) << delimiter;
      motif_file << tss_size << delimiter;
      motif_file << tss_preface_coding ->getSequenceAsString() << delimiter;
      motif_file << tss_coding->getSequenceAsString() << delimiter << '\n';

    }

    ExecEnv::log().info("Contig: {} has {} TSS blocks defined, assigned: {}", contig.first, vector.size(), assigned_count);

  }


}


void kgl::PromoterMotif::headerTFFMotif(std::ofstream& motif_file, const char delimiter) {

  motif_file << "GENE" << delimiter;
  motif_file << "Description" << delimiter;
  motif_file << "Symblic" << delimiter;
  motif_file << "AltSymbolic" << delimiter;
  motif_file << "Begin" << delimiter;
  motif_file << "End" << delimiter;
  motif_file << "Strand" << delimiter;
  motif_file << "5PrimeOffset" << delimiter;
  motif_file << "Offset" << delimiter;
  motif_file << "TFF" << delimiter;
  motif_file << "TFFBegin" << delimiter;
  motif_file << "TFFEnd" << delimiter;
  motif_file << "TFFStrand" << delimiter;
  motif_file << "TFFSize" << delimiter;
  motif_file << "TFFPreface" << delimiter;
  motif_file << "TFFSequence" << delimiter << '\n';

}