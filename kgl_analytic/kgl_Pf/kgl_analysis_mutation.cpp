//
// Created by kellerberrin on 5/1/21.
//

#include "kgl_analysis_mutation.h"
#include "kgl_analysis_gene_sequence.h"


namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kgl::MutationAnalysis::initializeAnalysis(const std::string& work_directory,
                                               const ActiveParameterList& named_parameters,
                                               std::shared_ptr<const GenomeCollection> reference_genomes) {

  ExecEnv::log().info("Default Analysis Id: {} initialized with work directory: {}", ident(), work_directory);
  for (auto const& [parameter_ident, parameter_map] : named_parameters.getMap()) {

    ExecEnv::log().info("Default Initialize Analysis Id: {}, initialized with parameter block: {}", ident(), parameter_ident);

  }

  if (not getParameters(named_parameters, reference_genomes, work_directory)) {

    ExecEnv::log().critical("Analysis Id: {}, problem parsing parameters, program ends. {}", ident());

  }

  return true;

}




bool kgl::MutationAnalysis::getParameters(const ActiveParameterList& named_parameters,
                                          const std::shared_ptr<const GenomeCollection>& genome_collection_ptr,
                                          const std::string& work_directory) {

  work_directory_ = work_directory;

  for (auto const& named_block : named_parameters.getMap()) {

    auto [block_name, block_vector] = named_block.second;

    if (block_vector.size() != 1) {

      ExecEnv::log().error("MutationAnalysis::getParameters; parameter block: {} vector size: {}, expected size = 1",
                           block_name, block_vector.size());
      return false;

    }

    ExecEnv::log().info("Analysis: {} parsing parameter block: {}", ident(), block_name);

    for (auto const& xml_vector : block_vector) {


      auto output_opt = xml_vector.getString(OUTPUT_FILE_);
      if (output_opt) {

        output_file_name_ = output_opt.value().front() + std::string(OUTPUT_FILE_EXT_);
        output_file_name_ = Utility::filePath(output_file_name_, work_directory);
        ExecEnv::log().info("Analysis: {} outputfile: {}", ident(), output_file_name_);

      } else {

        ExecEnv::log().error("MutationAnalysis::getParameters; bad value for parameter: {}", OUTPUT_FILE_);
        return false;

      }

      std::shared_ptr<GenomeCollection> genome_collection(std::make_shared<GenomeCollection>());
      genome_collection_ptr_ = genome_collection;

      // Get the vector of active genomes.
      auto genome_vector_opt = xml_vector.getString(GENOME_IDENT_, ParameterMap::ANY_SIZE);
      if (genome_vector_opt) {

        std::vector<GenomeId_t> genomes = genome_vector_opt.value();
        for (auto const& genome_id : genomes) {

          auto genome_opt = genome_collection_ptr->getOptionalGenome(genome_id);
          if (genome_opt) {

            ExecEnv::log().info("Analysis: {} add active genome: {}", ident(), genome_opt.value()->genomeId());
            if (not genome_collection->addGenome(genome_opt.value())) {

              ExecEnv::log().error("Analysis: {} active genome: {} cannot be added to active genome collection", ident(), genome_id);
              return false;

            }

          } else {

            ExecEnv::log().error("Analysis: {} active genome: {} not found in supplied genome collection", ident(), genome_id);
            return false;

          }

        } // for genomes.

        genome_collection_ptr_ = genome_collection;

      } else {

        ExecEnv::log().error("IMutationAnalysis::getParameters; bad value for parameter: {}", GENOME_IDENT_);
        return false;

      }

    }

  }

  return true;

}



// Perform the genetic analysis per iteration.
bool kgl::MutationAnalysis::writeOutput() {

  std::ofstream out_file(output_file_name_);

  if (not out_file.good()) {

    ExecEnv::log().error("MutationAnalysis::writeOutput; could not open file: {} for output", output_file_name_);
    return false;

  } else {

    ExecEnv::log().info("Analysis: {} writing to file: {} for output", ident(), output_file_name_);

  }

  writeHeader(out_file);

  for (auto const& [genome_id, genome_ptr] : genome_collection_ptr_->getMap()) {

    for (auto const& [contig_id, contig_ptr] : genome_ptr->getMap()) {

      for (auto const& [offset, gene_ptr] : contig_ptr->getGeneMap()) {

        std::shared_ptr<const CodingSequenceArray> sequence_array = GeneFeature::getCodingSequences(gene_ptr);

        out_file << genome_id << OUTPUT_DELIMITER_ << contig_id << OUTPUT_DELIMITER_
                 << gene_ptr->id() << OUTPUT_DELIMITER_ << offset << OUTPUT_DELIMITER_
                 << gene_ptr->sequence().begin() << OUTPUT_DELIMITER_
                 << gene_ptr->sequence().end() << OUTPUT_DELIMITER_
                 << gene_ptr->sequence().length() << OUTPUT_DELIMITER_
                 << gene_ptr->sequence().strandText() << OUTPUT_DELIMITER_
                 << sequence_array->size() << OUTPUT_DELIMITER_;

        auto const& attribute_map = gene_ptr->getAttributes().getMap();

        out_file << attribute_map.size();

        for (auto const& [attrib_key, attrib_value] : attribute_map) {

          out_file << OUTPUT_DELIMITER_ << attrib_key << ':' << attrib_value;

        }

        out_file << '\n';

      } // Gene.

    } // Contig.

  } // Genome

  return true;

}


void kgl::MutationAnalysis::writeHeader(std::ostream& out_file) {

  out_file << "Genome" << OUTPUT_DELIMITER_ << "Contig" << OUTPUT_DELIMITER_
           << "Gene" << OUTPUT_DELIMITER_ << "Offset" << OUTPUT_DELIMITER_
           << "Begin" << OUTPUT_DELIMITER_ << "End" << OUTPUT_DELIMITER_
           << "Length" << OUTPUT_DELIMITER_ << "Strand" << OUTPUT_DELIMITER_
           << "Exons" << OUTPUT_DELIMITER_ << "Attributes" << OUTPUT_DELIMITER_
           << "AttKey:AttVal" << '\n';

}


// Perform the genetic analysis per iteration.
bool kgl::MutationAnalysis::fileReadAnalysis(std::shared_ptr<const DataDB> data_ptr) {

  ExecEnv::log().info("File Read for Analysis Id: {} called with file: {}", ident(), data_ptr->fileId());

  auto file_characteristic = data_ptr->dataCharacteristic();

  if (file_characteristic.data_source == DataSourceEnum::Pf3kCOI) {

    pf3k_coi_ptr_ = std::dynamic_pointer_cast<const Pf3kCOIDB>(data_ptr);

    if (pf3k_coi_ptr_) {

      ExecEnv::log().info("Processed Pf3k complexity of infection file: {}", pf3k_coi_ptr_->fileId());

    } else {

      ExecEnv::log().critical("MutationAnalysis::fileReadAnalysis; Unable to cast to the Pf3K complexity of infection file, severe error.");

    }

  }


  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::MutationAnalysis::iterationAnalysis() {

  ExecEnv::log().info("Default Iteration Analysis called for Analysis Id: {}", ident());

  return true;

}

// All VCF data has been presented, finalize analysis and write results.
bool kgl::MutationAnalysis::finalizeAnalysis() {

  ExecEnv::log().info("Default Finalize Analysis called for Analysis Id: {}", ident());

  writeOutput();

  return true;

}



void kgl::MutationAnalysis::performRegion() {

  const GenomeId_t analysis_genome = "Pf3D7_47";

  std::string region_fasta_file = "malawi_fb_SRR609075";
  region_fasta_file += "_561666";
  region_fasta_file = Utility::filePath(region_fasta_file, work_directory_) + ".fasta";
  if (not GenomicSequence::mutateGenomeRegion("malawi_fb_SRR609075",
                                                   "Pf3D7_04_v3",
                                                   561666,
                                                   11753,
                                                   population_ptr_,
                                                   genome_collection_ptr_->getGenome(analysis_genome),
                                                   region_fasta_file)) {

    ExecEnv::log().error("PhylogeneticAnalysis::performRegion(); analysis fails");

  }

  region_fasta_file = "malawi_fb_SRR609075";
  region_fasta_file += "_462000";
  region_fasta_file = Utility::filePath(region_fasta_file, work_directory_) + ".fasta";
  if (not GenomicSequence::mutateGenomeRegion("malawi_fb_SRR609075",
                                                   "Pf3D7_01_v3",
                                                   462000,
                                                   2000,
                                                   population_ptr_,
                                                   genome_collection_ptr_->getGenome(analysis_genome),
                                                   region_fasta_file)) {

    ExecEnv::log().error("PhylogeneticAnalysis::performRegion(); analysis fails");

  }

  region_fasta_file = "malawi_fb_SRR609075";
  region_fasta_file += "_0";
  region_fasta_file = Utility::filePath(region_fasta_file, work_directory_) + ".fasta";
  if (not GenomicSequence::mutateGenomeRegion("malawi_fb_SRR609075",
                                                   "Pf3D7_04_v3",
                                                   0,
                                                   1200490,
                                                   population_ptr_,
                                                   genome_collection_ptr_->getGenome(analysis_genome),
                                                   region_fasta_file)) {

    ExecEnv::log().error("PhylogeneticAnalysis::performRegion(); analysis fails");

  }

  region_fasta_file = "malawi_fb_SRR609075";
  region_fasta_file += "_944950";
  region_fasta_file = Utility::filePath(region_fasta_file, work_directory_) + ".fasta";
  if (not GenomicSequence::mutateGenomeRegion("malawi_fb_SRR609075",
                                                   "Pf3D7_04_v3",
                                                   944950,
                                                   135,
                                                   population_ptr_,
                                                   genome_collection_ptr_->getGenome(analysis_genome),
                                                   region_fasta_file)) {

    ExecEnv::log().error("PhylogeneticAnalysis::performRegion(); analysis fails");

  }

  region_fasta_file = "malawi_fb_SRR609075";
  region_fasta_file += "_941875";
  region_fasta_file = Utility::filePath(region_fasta_file, work_directory_) + ".fasta";
  if (not GenomicSequence::mutateGenomeRegion("malawi_fb_SRR609075",
                                                   "Pf3D7_04_v3",
                                                   941875,
                                                   4293,
                                                   population_ptr_,
                                                   genome_collection_ptr_->getGenome(analysis_genome),
                                                   region_fasta_file)) {

    ExecEnv::log().error("PhylogeneticAnalysis::performRegion(); analysis fails");

  }

  region_fasta_file = "malawi_fb_SRR609075";
  region_fasta_file += "_569342";
  region_fasta_file = Utility::filePath(region_fasta_file, work_directory_) + ".fasta";
  if (not GenomicSequence::mutateGenomeRegion("malawi_fb_SRR609075",
                                                   "Pf3D7_04_v3",
                                                   569342,
                                                   7467,
                                                   population_ptr_,
                                                   genome_collection_ptr_->getGenome(analysis_genome),
                                                   region_fasta_file)) {

    ExecEnv::log().error("PhylogeneticAnalysis::performRegion(); analysis fails");

  }

  region_fasta_file = "malawi_fb_SRR609075";
  region_fasta_file += "_584668";
  region_fasta_file = Utility::filePath(region_fasta_file, work_directory_) + ".fasta";
  if (not GenomicSequence::mutateGenomeRegion("malawi_fb_SRR609075",
                                                   "Pf3D7_04_v3",
                                                   584668,
                                                   7280,
                                                   population_ptr_,
                                                   genome_collection_ptr_->getGenome(analysis_genome),
                                                   region_fasta_file)) {

    ExecEnv::log().error("PhylogeneticAnalysis::performRegion(); analysis fails");

  }


}

