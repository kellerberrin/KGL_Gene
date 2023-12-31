//
// Created by kellerberrin on 3/11/23.
//


#include "kga_analysis_lib_Pfgene.h"
#include "kgl_distance_tree_upgma.h"
#include "kgl_distance_tree_base.h"
#include "kgl_sequence_distance_impl.h"
#include "kgl_sequence_node.h"
#include "kel_utility.h"


namespace kga = kellerberrin::genome::analysis;
namespace kgl = kellerberrin::genome;



void kga::AnalysisGenePf::performGeneAnalysis(const std::shared_ptr<const GenomeReference>& genome_ref_ptr,
                                              const std::string& ident_work_directory_) {

  const std::string& genome_ref_id = genome_ref_ptr->genomeId();

  std::string intron_file_name = INTRON_ + genome_ref_id + INTRON_EXT_;
  intron_file_name = Utility::filePath(intron_file_name, ident_work_directory_);

  auto var_gene_vector = getGeneVector(genome_ref_ptr, PFEMP1_FAMILY_);
  auto ruf6_gene_vector = getGeneVector(genome_ref_ptr, RUF6_FAMILY_);
  auto rifin_gene_vector = getGeneVector(genome_ref_ptr, RIFIN_FAMILY_);
  auto stevor_gene_vector = getGeneVector(genome_ref_ptr, STEVOR_FAMILY_);
  auto surfin_gene_vector = getGeneVector(genome_ref_ptr, SURFIN_FAMILY_);
  auto ncRNA_gene_vector = getncRNAGeneVector(genome_ref_ptr, RUF6_FAMILY_);

  varIntron(var_gene_vector, intron_file_name);

  std::string newick_file_name = NEWICK_ + genome_ref_id + NEWICK_EXT_;
  newick_file_name = Utility::filePath(newick_file_name, ident_work_directory_);

  geneCodingUPGMA(var_gene_vector, newick_file_name);

  ExecEnv::log().info("performGeneAnalysis, Var Genes: {}, Rifin Genes: {}, Stevor: {}, Surfin: {}, RUF6: {}, ncRNA: {}",
                      var_gene_vector.size(),
                      rifin_gene_vector.size(),
                      stevor_gene_vector.size(),
                      rifin_gene_vector.size(),
                      ruf6_gene_vector.size(),
                      ncRNA_gene_vector.size());

  newick_file_name = NEWICK_ + std::string(RUF6_FAMILY_) + genome_ref_id + NEWICK_EXT_;
  newick_file_name = Utility::filePath(newick_file_name, ident_work_directory_);

  for (auto const& gene_ptr : ncRNA_gene_vector) {

    size_t radius{25000};
    ExecEnv::log().info("PerformGeneAnalysis, genome: {}, ncRNA gene: {}, genetype: {}, radius: {}, var: {}, rifin: {} surfin: {} surfin: {}"
    , genome_ref_id, gene_ptr->id(), gene_ptr->type(), radius
    , proximityGenes(radius, gene_ptr, var_gene_vector).size()
    , proximityGenes(radius, gene_ptr, rifin_gene_vector).size()
    , proximityGenes(radius, gene_ptr, surfin_gene_vector).size()
    , proximityGenes(radius, gene_ptr, stevor_gene_vector).size());

    ExecEnv::log().info("performGeneAnalysis, ncRNA feature: {}", gene_ptr->featureText());

  }

  geneCodingUPGMA(ruf6_gene_vector, newick_file_name);

}


kgl::GeneVector kga::AnalysisGenePf::getGeneVector(const std::shared_ptr<const GenomeReference>& genome_ptr,
                                                   const std::string& desc_uc_text) {

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


kgl::GeneVector kga::AnalysisGenePf::getncRNAGeneVector(const std::shared_ptr<const GenomeReference>& genome_ptr,
                                                        const std::string& desc_uc_text,
                                                        size_t max_size) {

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


kgl::GeneVector kga::AnalysisGenePf::proximityGenes(size_t radius,
                                                    const std::shared_ptr<const GeneFeature>& target_ptr,
                                                    const GeneVector& gene_vector) {

  GeneVector same_contig;
  for (auto const& gene_ptr : gene_vector) {

    if (gene_ptr->contig_ref_ptr()->contigId() == target_ptr->contig_ref_ptr()->contigId()) {

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



void kga::AnalysisGenePf::varIntron(const GeneVector& gene_vector,
                                    const std::string& intron_file_name) {

  // Open a file to receive csv delimited intron information.
  std::ofstream intron(intron_file_name);
  if (not intron.good()) {

    ExecEnv::log().critical("Unable to open Intron data file: {}", intron_file_name);

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
         << i_promoter.getStringView() << CSV_DELIMITER_
         << i_complement_promoter.getStringView() << CSV_DELIMITER_
         << i_5_promoter.getStringView() << CSV_DELIMITER_
         << "Size" << CSV_DELIMITER_
         << "ISequence" << '\n';

  for (auto const& gene_ptr : gene_vector) {

    auto coding_seq_ptr = GeneFeature::getTranscriptionSequences(gene_ptr);
    if (coding_seq_ptr->empty()) {

      ExecEnv::log().critical("Gene contains no coding sequence : gene: {}", gene_ptr->id());

    }

    std::shared_ptr<const TranscriptionSequence> transcript_ptr = coding_seq_ptr->getFirst();
    auto coding_dna_opt = gene_ptr->contig_ref_ptr()->codingSequence(transcript_ptr);
    if (coding_dna_opt) {

      DNA5SequenceCoding& coding_dna_sequence = coding_dna_opt.value();
      // Do we have a valid intron (VAR only)?

      auto intron_map = transcript_ptr->getIntronIntervals();
      const auto& contig_sequence = gene_ptr->contig_ref_ptr()->sequence();

      // Only add genes with valid coding sequences (no pseudo genes).
      auto sequence_validity = gene_ptr->contig_ref_ptr()->checkValidCodingSequence(coding_dna_sequence);
      if (TranscriptionSequence::checkValidProtein(sequence_validity)) {

        // Only 1 intron (var genes)
        if (intron_map.size() == 1) {

          auto first_intron_opt = contig_sequence.subSequence(*intron_map.begin());
          if (not first_intron_opt) {

            ExecEnv::log().error("Failed to generate intron interval: {} for gene: {}",
                                 intron_map.begin()->toString(), gene_ptr->id());
            return;

          }
          const auto& intron_sequence_linear = first_intron_opt.value();
          // Get the intron in the coding sense.
          StrandSense strand = transcript_ptr->strand();
          auto intron_sequence_coding = intron_sequence_linear.codingSequence(strand);
          // Down convert back to linear DNA5 (retains coding sense).
          auto intron_sequence = DNA5SequenceLinear::downConvertToLinear(intron_sequence_coding);

          std::stringstream pss;
          std::vector<ContigOffset_t> offset_vec = intron_sequence.findAll(i_promoter);
          for (auto offset : offset_vec) {

            pss << offset << ":";

          }

          std::stringstream css;
          offset_vec = intron_sequence.findAll(i_complement_promoter);
          for (auto offset : offset_vec) {

            css << offset << ":";

          }

          std::stringstream fss;
          offset_vec = intron_sequence.findAll(i_5_promoter);
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
                 << static_cast<char>(strand) << CSV_DELIMITER_
                 << pss.str()
                 << CSV_DELIMITER_
                 << css.str()
                 << CSV_DELIMITER_
                 << fss.str()
                 << CSV_DELIMITER_
                 << intron_sequence.length()
                 << CSV_DELIMITER_
                 << intron_sequence.getStringView() << '\n';

        } // 1 intron

      } // Valid Gene.

    } // Get coding sequence

  } // Is family member.

  ExecEnv::log().info("varIntron, processed genes: {}", gene_vector.size());

}



// Function to compare a family of reference genes (unmutated)
void kga::AnalysisGenePf::geneCodingUPGMA(const GeneVector& gene_vector,
                                          const std::string& newick_file_name) {

  CodingDistanceMetric coding_distance_metric{LevenshteinLocalCoding};
  TreeNodeVector distance_node_vector;

  for (auto const& gene_ptr : gene_vector) {

    auto coding_seq_ptr = GeneFeature::getTranscriptionSequences(gene_ptr);
    if (coding_seq_ptr->empty()) {

      ExecEnv::log().error("Gene contains no coding sequence : gene: {}", gene_ptr->id());
      continue;

    }

    for (auto& [transcript_id, transcript_ptr] : coding_seq_ptr->getMap()) {

      auto coding_dna_opt = gene_ptr->contig_ref_ptr()->codingSequence(transcript_ptr);
      if (coding_dna_opt) {

        DNA5SequenceCoding& coding_dna_sequence = coding_dna_opt.value();
        auto coding_sequence_ptr = std::make_shared<CodingSequenceNode>(std::move(coding_dna_sequence), transcript_id, coding_distance_metric);
        distance_node_vector.push_back(coding_sequence_ptr);

      } // Valid coding sequence

    }

  } // All Genes.

  DistanceTreeUPGMA upgma_distance;
  // Add the parentDistance vector.
  upgma_distance.addDistanceMap(distance_node_vector);
  // Calculate.
  auto root_node_vector = upgma_distance.calculateTree();
  // Report Results.
  DistanceTreeBase upgma_tree(root_node_vector);
  if (not upgma_tree.writeNewick(newick_file_name)) {

    ExecEnv::log().error("Unable to write UPGMA Newick file: {}", newick_file_name);

  }

}

