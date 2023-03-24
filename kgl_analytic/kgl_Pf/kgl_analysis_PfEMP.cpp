//
// Created by kellerberrin on 3/1/21.
//

#include "kgl_upgma.h"

#include "kgl_analysis_PfEMP.h"


namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kgl::PfEMPAnalysis::initializeAnalysis(const std::string& work_directory,
                                            const ActiveParameterList& named_parameters,
                                            const std::shared_ptr<const AnalysisResources>& resource_ptr) {

  ExecEnv::log().info("Default Analysis Id: {} initialized with work directory: {}", ident(), work_directory);
  for (auto const& [parameter_ident, parameter_map] : named_parameters.getMap()) {

    ExecEnv::log().info("Default Initialize Analysis Id: {}, initialized with parameter: {}", ident(), parameter_ident);

  }

  // Setup and clear the directories to hold analysis output.
  // The top level directory for this analysis type.
  ident_work_directory_ = work_directory + std::string("/") + ident();
  if (not Utility::createDirectory(ident_work_directory_)) {

    ExecEnv::log().critical("PfEMPAnalysis::initializeAnalysis, unable to create analysis results directory: {}",
                            ident_work_directory_);

  }

  for (auto const& genome_resource_ptr : resource_ptr->getResources(RuntimeResourceType::GENOME_DATABASE)) {

    auto genome_ptr = std::dynamic_pointer_cast<const GenomeReference>(genome_resource_ptr);
    ExecEnv::log().info("Initialize for Analysis Id: {} called with Reference Genome: {}", ident(), genome_ptr->genomeId());
    reference_genomes_->addGenome(genome_ptr);

  }

  if (not getParameters(named_parameters)) {

    ExecEnv::log().info("PfEMPAnalysis::initializeAnalysis; Analysis Id: {} problem parsing parameters", ident());
    return false;

  }

  performPFEMP1UPGMA();

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::PfEMPAnalysis::fileReadAnalysis(std::shared_ptr<const DataDB> base_file_ptr) {

  ExecEnv::log().info("File Read for Analysis Id: {} called for file: {}", ident(), base_file_ptr->fileId());

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::PfEMPAnalysis::iterationAnalysis() {

  ExecEnv::log().info("Default Iteration Analysis called for Analysis Id: {}", ident());

  return true;

}

// All VCF data has been presented, finalize analysis and write results.
bool kgl::PfEMPAnalysis::finalizeAnalysis() {

  ExecEnv::log().info("Default Finalize Analysis called for Analysis Id: {}", ident());

  return true;

}



void kgl::PfEMPAnalysis::performPFEMP1UPGMA() {

  
  VarGeneFamilyTree( newick_file_name_,
                     intron_file_name_,
                     reference_genomes_,
                     PFEMP1_FAMILY_);

}



bool kgl::PfEMPAnalysis::getParameters(const ActiveParameterList& named_parameters) {


  for (auto const& named_block : named_parameters.getMap()) {

    auto [block_name, block_vector] = named_block.second;

    if (block_vector.size() != 1) {

      ExecEnv::log().error("PfEMPAnalysis::getParameters; parameter block: {} vector size: {}, expected size = 1",
                           block_name, block_vector.size());
      return false;

    }

    ExecEnv::log().info("Analysis: {} parsing parameter block: {}", ident(), block_name);

    for (auto const& xml_vector : block_vector) {


      auto newick_opt = xml_vector.getString(NEWICK_FILE_);
      if (newick_opt) {

        newick_file_name_ = newick_opt.value().front();
        newick_file_name_ = Utility::filePath(newick_file_name_, ident_work_directory_);

      } else {

        ExecEnv::log().error("PfEMPAnalysis::getParameters; bad value for parameter: {}", NEWICK_FILE_);
        return false;

      }

      auto intron_opt = xml_vector.getString(INTRON_FILE_);
      if (intron_opt) {

        intron_file_name_ = intron_opt.value().front();
        intron_file_name_ = Utility::filePath(intron_file_name_, ident_work_directory_);

      } else {

        ExecEnv::log().error("PfEMPAnalysis::getParameters; bad value for parameter: {}", INTRON_FILE_);
        return false;

      }

    }

  }

  return true;

}



// Function (not variadic) to combine the UPGMAMatrix and UPGMADistanceNode to compare a family of reference genes (unmutated)
void kgl::PfEMPAnalysis::VarGeneFamilyTree( const std::string& newick_file,
                        const std::string& intron_file,
                        std::shared_ptr<const GenomeCollection> genome_collection_ptr,
                        const std::string& protein_family) {

  std::shared_ptr<const LevenshteinLocal> sequence_distance(std::make_shared<const LevenshteinLocal>());
  UPGMAMatrix upgma_distance;

  // Open a file to receive csv delimited intron information.
  const char delimiter = ',';
  std::ofstream intron(intron_file);

  if (not intron.good()) {

    ExecEnv::log().critical("ReferenceGeneDistance::getSequence(); Unable to open Intron data file: {}", intron_file);

  }

#define I_PROMOTER "TGTATGTG"
#define I_COMPLEMENT_PROMOTER "ACATACAC"
#define I_5_PROMOTER "TCATA"
  // Header.
  DNA5SequenceLinear i_promoter(StringDNA5(I_PROMOTER));
  DNA5SequenceLinear i_promoter_rev = DNA5SequenceLinear::downConvertToLinear(i_promoter.codingSequence(StrandSense::REVERSE));
  //  i_promoter_ptr = i_promoter_ptr_rev;

  DNA5SequenceLinear i_complement_promoter(StringDNA5(I_COMPLEMENT_PROMOTER));
  DNA5SequenceLinear i_complement_promoter_rev = DNA5SequenceLinear::downConvertToLinear(i_complement_promoter.codingSequence(StrandSense::REVERSE));
  //  i_complement_promoter_ptr = i_complement_promoter_ptr_rev;

  DNA5SequenceLinear i_5_promoter(StringDNA5(I_5_PROMOTER));
  DNA5SequenceLinear i_5_promoter_rev = DNA5SequenceLinear::downConvertToLinear(i_5_promoter.codingSequence(StrandSense::REVERSE));
  //  i_5_promoter_ptr = i_5_promoter_ptr_rev;

  intron << "Gene" << delimiter
         << "description" << delimiter
         << "IStart" << delimiter
         << "IEnd" << delimiter
         << "IStrand" << delimiter
         << i_promoter.getSequenceAsString() << delimiter
         << i_complement_promoter.getSequenceAsString() << delimiter
         << i_5_promoter.getSequenceAsString() << delimiter
         << "Size" << delimiter
         << "ISequence" << '\n';

  std::shared_ptr<PhyloNodeVector> node_vector_ptr(std::make_shared<PhyloNodeVector>());
  size_t gene_count{0};
  for (auto const& [genome_ident, genome_ptr] : genome_collection_ptr->getMap()) {

    for (auto const& [contig_ident, contig_ptr] : genome_ptr->getMap()) {

      ExecEnv::log().info("VarGeneFamilyTree; Contig: {} contains genes: {}", contig_ident, contig_ptr->getGeneMap().size());

      for (auto const& [gene_offset, gene_ptr] : contig_ptr->getGeneMap()) {

        auto description_vector = gene_ptr->getAttributes().getDescription();
        if (description_vector.empty()) {

          continue;

        }
        std::string description = description_vector.front();

        // Is this gene a member of the requested family.
        if (Utility::toupper(description).find(protein_family) != std::string::npos) {

          ++gene_count;

          const std::shared_ptr<const CodingSequenceArray> coding_seq_ptr = GeneFeature::getCodingSequences(gene_ptr);

          if (coding_seq_ptr->empty()) {

            ExecEnv::log().critical("ReferenceGeneDistance::getSequence(); Gene contains no coding sequence : gene: {}", gene_ptr->id());

          }

          std::shared_ptr<const CodingSequence> coding_sequence = coding_seq_ptr->getFirst();
          DNA5SequenceCoding coding_dna_sequence;

          if (contig_ptr->getDNA5SequenceCoding(coding_sequence, coding_dna_sequence)) {

            // Do we have a valid intron (VAR only)?
            std::vector<DNA5SequenceCoding> intron_sequence_array = contig_ptr->sequence().intronArraySequence(coding_sequence);
            std::shared_ptr<AminoGeneDistance> distance_ptr(std::make_shared<AminoGeneDistance>(sequence_distance,
                                                                                                genome_ptr,
                                                                                                gene_ptr,
                                                                                                protein_family));

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
            if (contig_ptr->verifyDNACodingSequence(coding_dna_sequence)) {

              // Only 1 intron (var genes)
              if (intron_sequence_array.size() == 1) {


                auto intron_sequence_iter = intron_sequence_array.begin();


                for (auto intron_seq : intron_offset_map) {

                  DNA5SequenceLinear intron_ptr = DNA5SequenceLinear::downConvertToLinear(*intron_sequence_iter);
                  DNA5SequenceLinear intron_seq_ptr_rev = DNA5SequenceLinear::downConvertToLinear(intron_ptr.codingSequence(StrandSense::REVERSE));
                  DNA5SequenceLinear intron_seq_ptr = DNA5SequenceLinear::strandedDownConvert(*intron_sequence_iter);
                  //                  intron_seq_ptr = intron_ptr;


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

                  intron << gene_ptr->id() << delimiter
                         << description << delimiter
                         << intron_seq.first << delimiter
                         << intron_seq.second << delimiter
                         << static_cast<char>(strand) << delimiter
                         << pss.str()
                         << delimiter
                         << css.str()
                         << delimiter
                         << fss.str()
                         << delimiter
                         << intron_seq_ptr.length()
                         << delimiter
                         << intron_seq_ptr.getSequenceAsString() << '\n';

                  ++intron_sequence_iter;

                }


#define MIN_INTRON_LENGTH 10
                if (intron_sequence_array.front().length() >= MIN_INTRON_LENGTH) {

                  std::shared_ptr<PhyloNode> phylo_node_ptr(std::make_shared<PhyloNode>(distance_ptr));
                  node_vector_ptr->push_back(phylo_node_ptr);

                } // min length

              } // 1 intron

            } // Valid Gene.

          } // Get coding sequence

        } // Is family member.

      } // All genes.

    } // All contigs.

  } // All Genomes

  ExecEnv::log().info("VarGeneFamilyTree; Gene family: {}, found genes: {}", protein_family, gene_count);
  // Calculate.
  upgma_distance.calculateTree(node_vector_ptr);
  // Report Results.
  upgma_distance.writeNewick(newick_file);

}

