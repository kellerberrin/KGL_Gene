//
// Created by kellerberrin on 6/04/23.
//


#include "kgl_upgma.h"
#include "kgl_analysis_PfEMP.h"


namespace kgl = kellerberrin::genome;



void kgl::PfEMPAnalysis::performPFEMP1UPGMA() {

  for (auto const& [genome_id, genome_ptr] : all_reference_genomes_ptr_->getMap()) {

    std::string intron_file_name = INTRON_ + genome_id + INTRON_EXT_;
    intron_file_name = Utility::filePath(intron_file_name, ident_work_directory_);

    auto var_gene_vector = getGeneVector(genome_ptr, PFEMP1_FAMILY_);
    auto ruf6_gene_vector = getGeneVector(genome_ptr, RUF6_FAMILY_);
    auto rifin_gene_vector = getGeneVector(genome_ptr, RIFIN_FAMILY_);
    auto stevor_gene_vector = getGeneVector(genome_ptr, STEVOR_FAMILY_);
    auto surfin_gene_vector = getGeneVector(genome_ptr, SURFIN_FAMILY_);
    auto ncRNA_gene_vector = getncRNAGeneVector(genome_ptr);

    varIntron(var_gene_vector, intron_file_name);

    std::string newick_file_name = NEWICK_ + genome_id + NEWICK_EXT_;
    newick_file_name = Utility::filePath(newick_file_name, ident_work_directory_);

    geneFamilyUPGMA(genome_ptr, var_gene_vector, newick_file_name, PFEMP1_FAMILY_);

    ExecEnv::log().info("PfEMPAnalysis::performPFEMP1UPGMA, Var Genes: {}, Rifin Genes: {}, Stevor: {}, Surfin: {}, RUF6: {}, ncRNA: {}",
                        var_gene_vector.size(), rifin_gene_vector.size(), stevor_gene_vector.size(),
                        rifin_gene_vector.size(), ruf6_gene_vector.size(), ncRNA_gene_vector.size());

    newick_file_name = NEWICK_ + std::string(RUF6_FAMILY_) + genome_id + NEWICK_EXT_;
    newick_file_name = Utility::filePath(newick_file_name, ident_work_directory_);

    for (auto const& gene_ptr : ncRNA_gene_vector) {

      size_t radius{25000};
      ExecEnv::log().info("PfEMPAnalysis::performPFEMP1UPGMA, genome: {}, ncRNA gene: {}, genetype: {}, radius: {}, var: {}, rifin: {} surfin: {} surfin: {}"
          , genome_id, gene_ptr->id(), gene_ptr->type(), radius
          , proximityGenes(radius, gene_ptr, var_gene_vector).size()
          , proximityGenes(radius, gene_ptr, rifin_gene_vector).size()
          , proximityGenes(radius, gene_ptr, surfin_gene_vector).size()
          , proximityGenes(radius, gene_ptr, stevor_gene_vector).size());

      ExecEnv::log().info("PfEMPAnalysis::performPFEMP1UPGMA, ncRNA feature: {}", gene_ptr->featureText());

    }

    //    geneFamilyUPGMA(genome_ptr, ruf6_gene_vector, newick_file_name, RUF6_FAMILY_);

  }


}


kgl::GeneVector kgl::PfEMPAnalysis::getGeneVector( const std::shared_ptr<const GenomeReference>& genome_ptr
    , const std::string& desc_uc_text) const {

  GeneVector gene_vector;

  for (auto const& [contig_ident, contig_ptr] : genome_ptr->getMap()) {

    for (auto const &[gene_offset, gene_ptr]: contig_ptr->getGeneMap()) {

      auto description_vector = gene_ptr->getAttributes().getDescription();
      for (auto const& description : description_vector) {

        if (Utility::toupper(description).find(desc_uc_text) != std::string::npos) {

          gene_vector.push_back(gene_ptr);

        }

      } // Gene

    } // Contig

  } // Genome

  return gene_vector;

}


kgl::GeneVector kgl::PfEMPAnalysis::getncRNAGeneVector( const std::shared_ptr<const GenomeReference>& genome_ptr,
                                                        const std::string& desc_uc_text,
                                                        size_t max_size) const {

  GeneVector gene_vector;
  size_t contig_count{0};

  for (auto const& [contig_ident, contig_ptr] : genome_ptr->getMap()) {

    ++contig_count;

    for (auto const &[gene_offset, gene_ptr]: contig_ptr->getGeneMap()) {

      auto description = gene_ptr->descriptionText();
      if (GeneFeature::ncRNACoding(gene_ptr)
          and  (gene_ptr->sequence().length() <= max_size or max_size == 0)
          and (Utility::toupper(description).find(desc_uc_text) != std::string::npos or desc_uc_text.empty())) {

        gene_vector.push_back(gene_ptr);

      }

    } // Gene

  } // Contig

  return gene_vector;

}



[[nodiscard]] kgl::GeneVector kgl::PfEMPAnalysis::proximityGenes(size_t radius,
                                                                 const std::shared_ptr<const GeneFeature>& target_ptr,
                                                                 const GeneVector& gene_vector) const {

  GeneVector same_contig;
  for (auto const& gene_ptr : gene_vector) {

    if (gene_ptr->contig()->contigId() == target_ptr->contig()->contigId()) {

      same_contig.push_back(gene_ptr);

    }
  }

  if (radius == 0) {

    return same_contig;

  }

  GeneVector nearby_genes;
  for (auto const& gene_ptr : same_contig) {


    if (target_ptr->sequence().distance(gene_ptr->sequence()) <= radius) {

      nearby_genes.push_back(gene_ptr);

    }

  }

  return nearby_genes;

}



void kgl::PfEMPAnalysis::varIntron( const GeneVector& gene_vector,
                                    const std::string& intron_file_name) const {

  std::shared_ptr<const LevenshteinLocal> sequence_distance(std::make_shared<const LevenshteinLocal>());

  // Open a file to receive csv delimited intron information.
  std::ofstream intron(intron_file_name);
  if (not intron.good()) {

    ExecEnv::log().critical("PfEMPAnalysis::varIntron; Unable to open Intron data file: {}", intron_file_name);

  }

  // Header.
  StringDNA5 i_prom(I_PROMOTER_);
  DNA5SequenceLinear i_promoter(std::move(i_prom));
  DNA5SequenceLinear i_promoter_rev = DNA5SequenceLinear::downConvertToLinear(i_promoter.codingSequence(StrandSense::REVERSE));
  //  i_promoter_ptr = i_promoter_ptr_rev;

  //  StringDNA5
  DNA5SequenceLinear i_complement_promoter(std::move(StringDNA5(I_COMPLEMENT_PROMOTER_)));
  DNA5SequenceLinear i_complement_promoter_rev = DNA5SequenceLinear::downConvertToLinear(i_complement_promoter.codingSequence(StrandSense::REVERSE));
  //  i_complement_promoter_ptr = i_complement_promoter_ptr_rev;

  DNA5SequenceLinear i_5_promoter(std::move(StringDNA5(I_5_PROMOTER_)));
  DNA5SequenceLinear i_5_promoter_rev = DNA5SequenceLinear::downConvertToLinear(i_5_promoter.codingSequence(StrandSense::REVERSE));
  //  i_5_promoter_ptr = i_5_promoter_ptr_rev;

  intron << "Gene" << CSV_DELIMITER_
         << "description" << CSV_DELIMITER_
         << "IStart" << CSV_DELIMITER_
         << "IEnd" << CSV_DELIMITER_
         << "IStrand" << CSV_DELIMITER_
         << i_promoter.getSequenceAsString() << CSV_DELIMITER_
         << i_complement_promoter.getSequenceAsString() << CSV_DELIMITER_
         << i_5_promoter.getSequenceAsString() << CSV_DELIMITER_
         << "Size" << CSV_DELIMITER_
         << "ISequence" << '\n';

  for (auto const& gene_ptr : gene_vector) {

    auto coding_seq_ptr = GeneFeature::getTranscriptionSequences(gene_ptr);

    if (coding_seq_ptr->empty()) {

      ExecEnv::log().critical("ReferenceGeneDistance::getSequence(); Gene contains no coding sequence : gene: {}", gene_ptr->id());

    }

    std::shared_ptr<const TranscriptionSequence> coding_sequence = coding_seq_ptr->getFirst();
    DNA5SequenceCoding coding_dna_sequence;

    if (gene_ptr->contig()->getDNA5SequenceCoding(coding_sequence, coding_dna_sequence)) {

      // Do we have a valid intron (VAR only)?
      std::vector<DNA5SequenceCoding> intron_sequence_array = gene_ptr->contig()->sequence().intronArraySequence(coding_sequence);
      StrandSense strand;
      IntronOffsetMap intron_offset_map;
      if (not SequenceOffset::intronOffsetAdapter(coding_sequence, strand, intron_offset_map)) {

        ExecEnv::log().error("VarGeneFamilyTree(), Contig: {}, Gene: {},  cannot generate INTRON map",
                             coding_sequence->contig()->contigId(), coding_sequence->getGene()->id());
        intron_offset_map.clear();
      }

      if (intron_offset_map.size() != intron_sequence_array.size()) {

        ExecEnv::log().error("UPGMAGeneFamilyTree, Intron map size: {} different from Inron sequence size: {}",
                             intron_offset_map.size(), intron_sequence_array.size());
        continue;

      }

      // Only add genes with valid coding sequences (no pseudo genes).
      if (gene_ptr->contig()->verifyDNACodingSequence(coding_dna_sequence)) {

        // Only 1 intron (var genes)
        if (intron_sequence_array.size() == 1) {

          auto& first_intron = intron_sequence_array.front();

          auto& intron_seq = *intron_offset_map.begin();

          DNA5SequenceLinear intron_ptr = DNA5SequenceLinear::downConvertToLinear(first_intron);
          DNA5SequenceLinear intron_seq_ptr_rev = DNA5SequenceLinear::downConvertToLinear(intron_ptr.codingSequence(StrandSense::REVERSE));
          DNA5SequenceLinear intron_seq_ptr = DNA5SequenceLinear::strandedDownConvert(first_intron);

          std::stringstream pss;
          std::vector<ContigOffset_t> offset_vec = intron_seq_ptr.findAll(i_promoter);
          for (auto offset : offset_vec) {

            pss << offset << ":";

          }

          std::stringstream css;
          offset_vec = intron_seq_ptr.findAll(i_complement_promoter);
          for (auto offset : offset_vec) {

            css << offset << ":";

          }

          std::stringstream fss;
          offset_vec = intron_seq_ptr.findAll(i_5_promoter);
          for (auto offset : offset_vec) {

            fss << offset << ":";

          }

          auto description_vector = gene_ptr->getAttributes().getDescription();
          std::string cat_description;
          for (auto const& desc : description_vector) {

            cat_description += desc;

          }

          intron << gene_ptr->id() << CSV_DELIMITER_
                 << cat_description << CSV_DELIMITER_
                 << intron_seq.first << CSV_DELIMITER_
                 << intron_seq.second << CSV_DELIMITER_
                 << static_cast<char>(strand) << CSV_DELIMITER_
                 << pss.str()
                 << CSV_DELIMITER_
                 << css.str()
                 << CSV_DELIMITER_
                 << fss.str()
                 << CSV_DELIMITER_
                 << intron_seq_ptr.length()
                 << CSV_DELIMITER_
                 << intron_seq_ptr.getSequenceAsString() << '\n';

        } // 1 intron

      } // Valid Gene.

    } // Get coding sequence

  } // Is family member.

  ExecEnv::log().info("PfEMPAnalysis::varIntron, processed genes: {}", gene_vector.size());

}



// Function (not variadic) to combine the UPGMAMatrix and UPGMADistanceNode to compare a family of reference genes (unmutated)
void kgl::PfEMPAnalysis::geneFamilyUPGMA( const std::shared_ptr<const GenomeReference>& genome_ptr,
                                          const GeneVector& gene_vector,
                                          const std::string& newick_file_name,
                                          const std::string& family_text) const {

  std::shared_ptr<const LevenshteinLocal> sequence_distance(std::make_shared<const LevenshteinLocal>());
  UPGMAMatrix upgma_distance;

  std::shared_ptr<PhyloNodeVector> node_vector_ptr(std::make_shared<PhyloNodeVector>());

  for (auto const& gene_ptr : gene_vector) {

    auto coding_seq_ptr = GeneFeature::getTranscriptionSequences(gene_ptr);

    if (coding_seq_ptr->empty()) {

      ExecEnv::log().critical("ReferenceGeneDistance::getSequence(); Gene contains no coding sequence : gene: {}", gene_ptr->id());

    }

    std::shared_ptr<const TranscriptionSequence> coding_sequence = coding_seq_ptr->getFirst();
    DNA5SequenceCoding coding_dna_sequence;

    if (gene_ptr->contig()->getDNA5SequenceCoding(coding_sequence, coding_dna_sequence)) {

      std::shared_ptr<AminoGeneDistance> distance_ptr(std::make_shared<AminoGeneDistance>(sequence_distance,
                                                                                          genome_ptr,
                                                                                          gene_ptr,
                                                                                          family_text));

      if (coding_dna_sequence.length() >= MIN_SEQUENCE_LENGTH_) {

        std::shared_ptr<PhyloNode> phylo_node_ptr(std::make_shared<PhyloNode>(distance_ptr));
        node_vector_ptr->push_back(phylo_node_ptr);

      } // min length

    } // 1 intron

  } // Valid Gene.

  // Calculate.
  upgma_distance.calculateTree(node_vector_ptr);
  // Report Results.
  if (not upgma_distance.writeNewick(newick_file_name)) {

    ExecEnv::log().error("VarGeneFamilyTree; Unable to write UPGMA Newick file: {}", newick_file_name);

  }

}


