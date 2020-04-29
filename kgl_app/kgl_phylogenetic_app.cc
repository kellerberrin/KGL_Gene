//
// Created by kellerberrin on 10/11/17.
//

#include "kgl_gff_fasta.h"
#include "kgl_variant_phase.h"
#include "kgl_phylogenetic_app.h"
#include "kgl_phylogenetic_app_analysis.h"
#include "kgl_distribution_analysis.h"
#include "kgl_filter.h"
#include "kgl_properties.h"
#include "kgl_variant_factory_vcf.h"


namespace kgl = kellerberrin::genome;



void kgl::PhylogeneticExecEnv::executeApp() {

  // Command line arguments
  const Phylogenetic& args = getArgs();
  // XML program run options
  const RuntimeProperties& runtime_options = getRuntimeOptions();

  // Create a collection of all active organism genomes.
  std::shared_ptr<kgl::GenomeCollection> genome_collection = GenomeCollection::createGenomeCollection(runtime_options);

  // Create a phased population object.
  std::shared_ptr<PhasedPopulation> population_ptr(std::make_shared<PhasedPopulation>("Organism"));

  // Create an unphased population object.
  std::shared_ptr<UnphasedPopulation> unphased_population_ptr(std::make_shared<UnphasedPopulation>("ParsedVCF"));

  // Write Filtered Unphased Heterozygous Statistics
  HeterozygousStatistics heterozygous_statistics;

  std::vector<VCFFileInfo> vcf_list = runtime_options.getVCFFileVector();
    // For all VCF files, read in the variants.
  for (const auto &vcf_file : vcf_list) {

    // Get VCF reference genome.
    std::shared_ptr<const GenomeDatabase> reference_genome_ptr = genome_collection->getGenome(
    vcf_file.referenceGenome());

    // Filter and process Gatk variants.
    if (vcf_file.parserType() == VCFParserEnum::GatkMultiGenome) {

      // Read variants.
      std::shared_ptr<UnphasedPopulation> parsed_variants = VcfFactory::gatkMultiGenomeVCFVariants(
      reference_genome_ptr, vcf_file.fileName());

      // Basic statistics to output
      // unphased_population_ptr->popStatistics();
      // Filter unphased variants for minimum read statistics.
      std::shared_ptr<UnphasedPopulation> filtered_unphased_ptr = parsed_variants->filterVariants(
      AndFilter(DPCountFilter(30), RefAltCountFilter(30)));

      // Process Filtered Unphased Heterozygous Statistics
      if (not heterozygous_statistics.heterozygousStatistics(filtered_unphased_ptr)) {

        ExecEnv::log().error(
        "PhylogeneticExecEnv::executeApp(), Cannot generate heterozygous statistics for VCF file: {}",
        vcf_file.fileName());

      }

      // If the mixture file is defined and exists then read it and generate a population of clonal genomes.
      std::string mixture_file;
      if (runtime_options.getMixtureFile(mixture_file)) {

        // The mixture file (Pf3k only) indicates the Complexity Of Infection (COI) of VCF samples, this function includes only clonal infections.
        std::shared_ptr<UnphasedPopulation> clonal_unphased = GenomePhasing::filterClonal(mixture_file, filtered_unphased_ptr);

      }

      // Get the VCF ploidy (need not be the organism ploidy).
      size_t ploidy = vcf_file.ploidy();
      // Phase the homozygous and heterozygous variants into a haploid population.
      GenomePhasing::haploidPhasing(ploidy, filtered_unphased_ptr, reference_genome_ptr, population_ptr);


      std::pair<size_t, size_t> valid_count = filtered_unphased_ptr->validate(reference_genome_ptr);
      ExecEnv::log().info("Population: {}, Total Variants: {}, Validated Variants: {}", filtered_unphased_ptr->populationId(), valid_count.first, valid_count.second);

      unphased_population_ptr->mergePopulation(filtered_unphased_ptr);

    } else if (vcf_file.parserType() == VCFParserEnum::GRChNoGenome) {

      ContigAliasMap contig_alias_map = runtime_options.getContigAlias();
      // Read variants.
      std::shared_ptr<UnphasedGenome> parsed_variants = VcfFactory::GRChNoGenomeVCFVariants(reference_genome_ptr,
                                                                                            vcf_file.fileName(),
                                                                                            contig_alias_map);

      std::pair<size_t, size_t> valid_count = parsed_variants->validate(reference_genome_ptr);
      ExecEnv::log().info("Genome: {}, Total Variants: {}, Validated Variants: {}", parsed_variants->genomeId(), valid_count.first, valid_count.second);

    }

  }


  // Write the unphased hetero/homozygous statistics.
  std::string heterozygous_file = Utility::filePath("HeterozygousAnalysis", args.workDirectory) + ".csv";
  if (not heterozygous_statistics.writeHeterozygousStatistics(heterozygous_file, ',')) {

    ExecEnv::log().error("PhylogeneticExecEnv::executeApp(), Cannot write heterozygous statistics to file: {}", heterozygous_file);

  }

  // Analyze the data.
  kgl::PhylogeneticAnalysis(runtime_options, genome_collection, unphased_population_ptr, population_ptr).performAnalysis(args.analysisType);

}

