//
// Created by kellerberrin on 3/1/21.
//

#include "kgl_analysis_PfEMP.h"
#include "kgl_variant_filter.h"
#include "kgl_variant_filter_db.h"
#include "kgl_variant_filter_info.h"


namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kgl::PfEMPAnalysis::initializeAnalysis(const std::string& work_directory,
                                            const ActiveParameterList& named_parameters,
                                            const std::shared_ptr<const AnalysisResources>& resource_ptr) {

  ExecEnv::log().info("Default Analysis Id: {} initialized with work directory: {}", ident(), work_directory);
  for (auto const& [parameter_ident, parameter_map] : named_parameters.getMap()) {

    ExecEnv::log().info("Default Initialize Analysis Id: {}, initialized with parameter: {}", ident(), parameter_ident);

  }

  Pf7_sample_ptr_ = resource_ptr->getSingleResource<const Pf7SampleResource>(ResourceProperties::PF7SAMPLE_RESOURCE_ID_);
  Pf7_fws_ptr_ = resource_ptr->getSingleResource<const Pf7FwsResource>(ResourceProperties::PF7FWS_RESOURCE_ID_);
  Pf7_distance_ptr_ = resource_ptr->getSingleResource<const Pf7DistanceResource>(ResourceProperties::PF7DISTANCE_RESOURCE_ID_);

  // Setup and clear the directories to hold analysis output.
  // The top level directory for this analysis type.
  ident_work_directory_ = work_directory + std::string("/") + ident();
  if (not Utility::createDirectory(ident_work_directory_)) {

    ExecEnv::log().critical("PfEMPAnalysis::initializeAnalysis, unable to create analysis results directory: {}",
                            ident_work_directory_);

  }

  auto reference_genomes_ptr = std::make_shared<GenomeCollection>();
  for (auto const& genome_resource_ptr : resource_ptr->getResources(ResourceProperties::GENOME_RESOURCE_ID_)) {

    auto genome_ptr = std::dynamic_pointer_cast<const GenomeReference>(genome_resource_ptr);
    ExecEnv::log().info("Initialize for Analysis Id: {} called with Reference Genome: {}", ident(), genome_ptr->genomeId());
    reference_genomes_ptr->addGenome(genome_ptr);

  }
  all_reference_genomes_ptr_ = reference_genomes_ptr;  // Assign to a pointer to const.

  auto pf3d7_opt = all_reference_genomes_ptr_->getOptionalGenome(PF3D7_IDENT_);
  if (not pf3d7_opt) {

    ExecEnv::log().critical("PfEMPAnalysis::initializeAnalysis; Reference Genome: {} required for analysis - not supplied", PF3D7_IDENT_);

  }
  genome_3D7_ptr_ = pf3d7_opt.value();

//  performPFEMP1UPGMA();

  // Select the genes we are interested in analyzing for genome variants.
  antigenic_gene_map_.setGeneVector(getAntiGenicGenes(genome_3D7_ptr_));
  all_gene_map_.setGeneVector(getAllGenes(genome_3D7_ptr_));
  translation_gene_map_.setGeneVector(getTranslationGenes(genome_3D7_ptr_));

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::PfEMPAnalysis::fileReadAnalysis(std::shared_ptr<const DataDB> base_data_ptr) {

  ExecEnv::log().info("File Read for Analysis Id: {} called for file: {}", ident(), base_data_ptr->fileId());

  // Superclass the population_ptr
  std::shared_ptr<const PopulationDB> population_ptr = std::dynamic_pointer_cast<const PopulationDB>(base_data_ptr);

  if (not population_ptr) {

    ExecEnv::log().error("Analysis: {}, expected a Population in file: {}", ident(), base_data_ptr->fileId());
    return false;

  }

  auto filtered_population_ptr = qualityFilter(population_ptr);

  ExecEnv::log().info("Population Returned Filtered Size Genome count: {}, Variant Count: {}",
                      filtered_population_ptr->getMap().size(),
                      filtered_population_ptr->variantCount());

//  checkDistanceMatrix(population_ptr, filtered_population_ptr);

  // Analyze gene variant info
  translation_gene_map_.getGeneVariants(filtered_population_ptr);
  antigenic_gene_map_.getGeneVariants(filtered_population_ptr);
  all_gene_map_.getGeneVariants(filtered_population_ptr);


  // Analyze for Homozygous and overlapping variants.
  hetero_homo_zygous_.analyzeVariantPopulation(filtered_population_ptr);

  return true;

}

// Quality filter the variants using read depth, VQSLOD and other statistics
std::shared_ptr<kgl::PopulationDB> kgl::PfEMPAnalysis::qualityFilter(const std::shared_ptr<const PopulationDB>& unfiltered_population_ptr) {


  size_t unfiltered_count = unfiltered_population_ptr->variantCount();
  ExecEnv::log().info("Unfiltered Population: {}, Genome count: {}, Variant Count: {}",
                      unfiltered_population_ptr->populationId(),
                      unfiltered_population_ptr->getMap().size(),
                      unfiltered_count);

  // Shallow filter only.
  auto filtered_qc_population_ptr = Pf7_sample_ptr_->filterPassQCGenomes(unfiltered_population_ptr);

  size_t qc_count = filtered_qc_population_ptr->variantCount();
  double qc_filtered = 100.0 * (static_cast<double>(unfiltered_count - qc_count) / static_cast<double>(unfiltered_count));
  ExecEnv::log().info("Filter P7 QC Pass Genome count: {}, Variant Count: {}, Sample Data Count: {}, %filtered: {}",
                      filtered_qc_population_ptr->getMap().size(),
                      qc_count,
                      Pf7_sample_ptr_->getMap().size(),
                      qc_filtered);

  // Shallow filter only.
  auto monoclonal_population_ptr = Pf7_fws_ptr_->filterFWS(FwsFilterType::GREATER_EQUAL, MONOCLONAL_FWS_THRESHOLD, filtered_qc_population_ptr);

  size_t monoclonal_count = monoclonal_population_ptr->variantCount();
  double mono_filtered = 100.0 * (static_cast<double>(qc_count - monoclonal_count) / static_cast<double>(qc_count));
  ExecEnv::log().info("Filter MonoClonal FWS: {}, Genome Count: {},Variant Count: {}, %filtered: {}",
                      MONOCLONAL_FWS_THRESHOLD,
                      monoclonal_population_ptr->getMap().size(),
                      monoclonal_count,
                      mono_filtered);

  // Filter on VQSLOD if it exists, if field not present then the variant is not filtered.
  auto vqslod_filter_lambda = [vsqlod=VQSLOD_LEVEL_](double compare) ->bool { return compare >= vsqlod; };
  auto vqslod_filter = InfoFilter<double, true>(VQSLOD_FIELD_, vqslod_filter_lambda);
  monoclonal_population_ptr->selfFilter(vqslod_filter);

  size_t vqslod_count = monoclonal_population_ptr->variantCount();
  double vqslod_filtered = 100.0 * (static_cast<double>(monoclonal_count - vqslod_count) / static_cast<double>(monoclonal_count));
  ExecEnv::log().info("Filter >= VQSLOD level {}, Variant Count: {}, %filtered: {}", VQSLOD_LEVEL_, vqslod_count, vqslod_filtered);

  // Filter on QD if it exists, if field not present then the variant is not filtered.
  auto qd_filter_lambda = [qd=QD_LEVEL_](double compare) ->bool { return compare >= qd; };
  auto qd_filter = InfoFilter<double, true>(QD_FIELD_, qd_filter_lambda);
  monoclonal_population_ptr->selfFilter(qd_filter);

  size_t qd_count = monoclonal_population_ptr->variantCount();
  double qd_filtered = 100.0 * (static_cast<double>(vqslod_count - qd_count) / static_cast<double>(vqslod_count));
  ExecEnv::log().info("Filter >= QD level {}, Variant Count: {}, %filtered: {}", QD_LEVEL_, qd_count, qd_filtered);

  // Filter on MQ if it exists, if field not present then the variant is not filtered.
  auto mq_filter_lambda = [mq=MQ_LEVEL_](double compare) ->bool { return compare >= mq; };
  auto mq_filter = InfoFilter<double, true>(MQ_FIELD_, mq_filter_lambda);
  monoclonal_population_ptr->selfFilter(mq_filter);

  size_t mq_count = monoclonal_population_ptr->variantCount();
  double mq_filtered = 100.0 * (static_cast<double>(qd_count - mq_count) / static_cast<double>(qd_count));
  ExecEnv::log().info("Filter >= MQ level {}, Variant Count: {}, %filtered: {}", MQ_LEVEL_, mq_count, mq_filtered);

  // Filter on SOR if it exists, if field not present then the variant is not filtered.
  auto sor_filter_lambda = [sor=SOR_LEVEL_](double compare) ->bool { return compare <= sor; };
  auto sor_filter = InfoFilter<double, true>(SOR_FIELD_, sor_filter_lambda);
  monoclonal_population_ptr->selfFilter(sor_filter);

  size_t sor_count = monoclonal_population_ptr->variantCount();
  double sor_filtered = 100.0 * (static_cast<double>(mq_count - sor_count) / static_cast<double>(mq_count));
  ExecEnv::log().info("Filter <= SOR level {}, Variant Count: {}, %filtered: {}", SOR_LEVEL_, sor_count, sor_filtered);

  // Filter on MQRankSum if it exists, if field not present then the variant is not filtered.
  auto mq_rank_sum_filter_lambda = [mqrs=MQRankSum_LEVEL_](double compare) ->bool { return compare >= mqrs; };
  auto mq_rank_sum_filter = InfoFilter<double, true>(MQRankSum_FIELD_, mq_rank_sum_filter_lambda);
  monoclonal_population_ptr->selfFilter(mq_rank_sum_filter);

  size_t mqrs_count = monoclonal_population_ptr->variantCount();
  double mqrs_filtered = 100.0 * (static_cast<double>(sor_count - mqrs_count) / static_cast<double>(sor_count));
  ExecEnv::log().info("Filter >= MQRankSum level {}, Variant Count: {}, %filtered: {}", MQRankSum_LEVEL_, mqrs_count, mqrs_filtered);

  // Filter on MQRankSum if it exists, if field not present then the variant is not filtered.
  auto read_sum_filter_lambda = [rprs=ReadPosRankSum_LEVEL_](double compare) ->bool { return compare >= rprs; };
  auto read_sum_filter = InfoFilter<double, true>(ReadPosRankSum_FIELD_, read_sum_filter_lambda);
  monoclonal_population_ptr->selfFilter(read_sum_filter);

  size_t rprs_count = monoclonal_population_ptr->variantCount();
  double rprs_filtered = 100.0 * (static_cast<double>(mqrs_count - rprs_count) / static_cast<double>(mqrs_count));
  ExecEnv::log().info("Filter >= ReadPosRankSum {}, Variant Count: {}, %filtered: {}", ReadPosRankSum_LEVEL_, rprs_count, rprs_filtered);

  // Filtered for variant quality. Read depth >= 10.
  monoclonal_population_ptr->selfFilter(DPCountFilter(MINIMUM_READ_DEPTH_));

  size_t depth_count = monoclonal_population_ptr->variantCount();
  double depth_filtered = 100.0 * (static_cast<double>(rprs_count - depth_count) / static_cast<double>(rprs_count));
  ExecEnv::log().info("Filter >= Read Depth: {} Genome count: {}, Variant Count: {}, %filtered: {}",
                      MINIMUM_READ_DEPTH_,
                      monoclonal_population_ptr->getMap().size(),
                      depth_count,
                      depth_filtered);

  // Filter for unique variants at each offset.
  //  filtered_monoclonal_ptr->selfFilter(UniqueUnphasedFilter());

  // We need to do a deep copy of the filtered population here since the pass QC and FWS P7 filters only do a shallow copy.
  // And when the resultant population pointers go out of scope they will take the shared population structure with them.
  auto filtered_population_ptr = monoclonal_population_ptr->deepCopy();
  // Filtered population should contain all contigs for all genomes.
  filtered_population_ptr->squareContigs();

  ExecEnv::log().info("Population Final Filtered Size Genome count: {}, Variant Count: {}",
                      filtered_population_ptr->getMap().size(),
                      filtered_population_ptr->variantCount());

  return filtered_population_ptr;

}

// Perform the genetic analysis per iteration.
bool kgl::PfEMPAnalysis::iterationAnalysis() {

  ExecEnv::log().info("Default Iteration Analysis called for Analysis Id: {}", ident());

  return true;

}

// All VCF data has been presented, finalize analysis and write results.
bool kgl::PfEMPAnalysis::finalizeAnalysis() {

  ExecEnv::log().info("Default Finalize Analysis called for Analysis Id: {}", ident());

  std::string variant_file_name = std::string(VARIANT_COUNT_) + "Translation" + std::string(VARIANT_COUNT_EXT_);
  variant_file_name = Utility::filePath(variant_file_name, ident_work_directory_);
  translation_gene_map_.writeGeneResults(variant_file_name);

  variant_file_name = std::string(VARIANT_COUNT_) + "AntiGenic" + std::string(VARIANT_COUNT_EXT_);
  variant_file_name = Utility::filePath(variant_file_name, ident_work_directory_);
  antigenic_gene_map_.writeGeneResults(variant_file_name);

  variant_file_name = std::string(VARIANT_COUNT_) + "All" + std::string(VARIANT_COUNT_EXT_);
  variant_file_name = Utility::filePath(variant_file_name, ident_work_directory_);
  all_gene_map_.writeGeneResults(variant_file_name);

  variant_file_name = std::string("VariantStatistics") + std::string(VARIANT_COUNT_EXT_);
  variant_file_name = Utility::filePath(variant_file_name, ident_work_directory_);
  hetero_homo_zygous_.write_results(variant_file_name);

  return true;

}



kgl::GeneVector kgl::PfEMPAnalysis::getAntiGenicGenes(const std::shared_ptr<const GenomeReference>& genome_ptr) {

  // Get the gene families of interest.
  auto var_gene_vector = getGeneVector(genome_ptr, PFEMP1_FAMILY_);
  auto ruf6_gene_vector = getGeneVector(genome_ptr, RUF6_FAMILY_);
  var_gene_vector.insert(var_gene_vector.end(), ruf6_gene_vector.begin(), ruf6_gene_vector.end() );
  auto rifin_gene_vector = getGeneVector(genome_ptr, RIFIN_FAMILY_);
  var_gene_vector.insert(var_gene_vector.end(), rifin_gene_vector.begin(), rifin_gene_vector.end() );
  auto stevor_gene_vector = getGeneVector(genome_ptr, STEVOR_FAMILY_);
  var_gene_vector.insert(var_gene_vector.end(), stevor_gene_vector.begin(), stevor_gene_vector.end() );
  auto surfin_gene_vector = getGeneVector(genome_ptr, SURFIN_FAMILY_);
  var_gene_vector.insert(var_gene_vector.end(), surfin_gene_vector.begin(), surfin_gene_vector.end() );
  //  auto ncRNA_gene_vector = getncRNAGeneVector(genome_ptr);
  //  var_gene_vector.insert( var_gene_vector.end(), ncRNA_gene_vector.begin(), ncRNA_gene_vector.end() );

  return var_gene_vector;

}

kgl::GeneVector kgl::PfEMPAnalysis::getTranslationGenes(const std::shared_ptr<const GenomeReference>& genome_ptr) {

  // Get the gene families of interest.
  auto trna_gene_vector = getncRNAGeneVector(genome_ptr, TRNA_FAMILY_);
  auto ribosome_gene_vector = getGeneVector(genome_ptr, RIBOSOME_FAMILY_);
  trna_gene_vector.insert(trna_gene_vector.end(), ribosome_gene_vector.begin(), ribosome_gene_vector.end() );

  return trna_gene_vector;

}


kgl::GeneVector kgl::PfEMPAnalysis::getAllGenes(const std::shared_ptr<const GenomeReference>& genome_ptr) {

  GeneVector all_genes;
  for (auto const& [contig_id, contig_ptr] : genome_ptr->getMap()) {

    for (auto const& [gene_id, gene_ptr] : contig_ptr->getGeneMap()) {

      all_genes.push_back(gene_ptr);

    }

  }

  return all_genes;

}

